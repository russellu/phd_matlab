clear all ; close all ; 
fmrivols = {'blur_reg_bp_mc_retino_gamma_01.nii.gz','blur_reg_bp_mc_retino_gamma_02.nii.gz'} ; 
stimtrigs = {'S  1','S  2','S  3'} ; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie'} ;
subcomps = {[28,24,20,32],[35,39],[22,23,52],[25,44],[27,26,29],[19,26,21],[44,21]} ;

for subby=1:7 ;
    clear minvals mininds corrvols diffvols taskvols basevols allts
    name = subs{subby} ;
    stimtrigcounts = [1,1,1] ;
    for eset = 1:length(fmrivols) 
    cd(['c:/shared/badger_eeg/',name,'/']) ; ls 
    setnames = dir('preproc*gamma*set') ; 
    eegsets = {setnames(1).name,setnames(2).name} ; 
    eeg = pop_loadset(eegsets{eset}) ; 
    cd(['C:\shared\badger_mri\',name,'\nii']) ; 
    bold = load_untouch_nii(fmrivols{eset}) ; 
    eegevents = {eeg.urevent.type} ; 
    eeglats = {eeg.urevent.latency} ;  
    volonsets = find(strcmp('R128',eegevents)) ; 
    volonsetlats = cell2mat(eeglats(volonsets)) ; 
    basetime = 2 ; baseTR = ceil(basetime/bold.hdr.dime.pixdim(5)) ; 
    tasktime = 12 ; taskTR = ceil(tasktime/bold.hdr.dime.pixdim(5)) ; 
    hrf = spm_hrf(bold.hdr.dime.pixdim(5)) ; 
    task = zeros(1,taskTR) ; task(1:round(5/bold.hdr.dime.pixdim(5))) = 1 ; task = conv(task,hrf) ; task = task(1:taskTR) ; 
    peakinds = find(task==max(task))-2:find(task==max(task))+2 ; 

    boldim = bold.img ; 
    for st=1:3
        stimonsets = find(strcmp(stimtrigs{st},eegevents)) ; 
        trigonsetlats = cell2mat(eeglats(stimonsets)) ; 
        for i=1:length(trigonsetlats) ; 
            absdiffs = abs(volonsetlats-trigonsetlats(i)) ; 
            rawdiffs = volonsetlats - trigonsetlats(i) ;  
            minind = find(absdiffs==min(absdiffs),1) ; 
            minvals(i) = rawdiffs(minind) ;
            mininds(i) = minind ; 
        end
        for i=1:length(mininds) ; 
           taskvoli = boldim(:,:,:,mininds(i):mininds(i)+taskTR-1) ;  
           corrvols(st,:,:,:,stimtrigcounts(st)) = voxcorr(taskvoli,task) ; 
           diffvols(st,:,:,:,stimtrigcounts(st)) = (squeeze(mean(taskvoli(:,:,:,peakinds),4)-mean(taskvoli(:,:,:,1:3),4)))./std(taskvoli,0,4) ; 
           taskvols(st,:,:,:,stimtrigcounts(st)) = squeeze(mean(taskvoli(:,:,:,peakinds),4)) ;
           basevols(st,:,:,:,stimtrigcounts(st)) = squeeze(mean(taskvoli(:,:,:,1:3),4)) ; 
           allts(st,stimtrigcounts(st),:,:,:,:) = taskvoli ;
           %mixvals(st,:,stimtrigcounts(st)) = squeeze(mixes(eset,mininds(i):mininds(i)+taskTR-1)) ;
           stimtrigcounts(st) = stimtrigcounts(st) + 1 ; 
        end
    end
    end
    restask = reshape(taskvols,[numel(taskvols(:,:,:,:,1)),32]) ; 
    resbase = reshape(basevols,[numel(basevols(:,:,:,:,1)),32]) ; 
    [h,p,ci,stats] = ttest(restask',resbase') ; 
    tstats = stats.tstat ; 
    tres = reshape(tstats,[3,size(corrvols,2),size(corrvols,3),size(corrvols,4)]) ; 
    pres = reshape(p,[3,size(corrvols,2),size(corrvols,3),size(corrvols,4)]) ;
    save('peakinds','peakinds') ; save('task','task') ; 
    save('allts','allts','-v7.3') ; 

    %%% get the subtracted power values

    comps = subcomps{subby}(1:2) ; allc = zeros(1,64) ; allc(comps) = 1 ; % first two components for back-projection
    elecs = [15,51,7,37,19,38,8,52,16,59,45,31,46,60,9,20,10] ; es = zeros(64,1) ; es(elecs) = 1 ; 
    clear baselow tasklow allepochs
    for eset=1:length(eegsets) 
        cd(['c:/shared/badger_eeg/',name,'/']) ; 
        eeg = pop_loadset(eegsets{eset}) ; 
        eeg = pop_subcomp(eeg,find(allc==0)) ; 
        fcount = 1 ; 
        for f=1:120
        eeglow = eeg ; eeglow.data = eegfiltfft(eeglow.data,eeglow.srate,f-2,f+2) ; 
        for st=1:3
            ep = pop_epoch(eeglow,{stimtrigs{st}},[-6,12]) ; 
            lowdata = ep.data.^2 ; 
    %        baselow(st,fcount,:,(eset*16-16)+1:eset*16) = squeeze(mean(lowdata(:,ep.times<0,:),2)) ; 
    %        tasklow(st,fcount,:,(eset*16-16)+1:eset*16) = squeeze(mean(lowdata(:,(ep.times>0) & ep.times<5000,:),2)) ; 
            allepochs(st,fcount,:,:,(eset*16-16)+1:eset*16) = squeeze(ep.data(elecs,:,:)) ; 
        end
        fcount = fcount + 1 ; 
        end
    end
    %lowdiff = tasklow-baselow ; 
    save('allepochs','allepochs','-v7.3') ; 

end

%{

% get posterior ROI
tthresh=7 ; 
roi = ((squeeze(mean(tres,1))>tthresh)) ; roivox = find(roi==1) ; 
for i=1:size(taskvols,1)
    for j=1:size(taskvols,5)
        boldij = squeeze(diffvols(i,:,:,:,j)) ; 
        boldtrials(i,j) = squeeze(mean(boldij(roivox))) ; 
    end
end
% correlate and plot
clear corrs ; 
for i=1:3
    for j=1:size(lowdiff,2) 
        corrs(i,j) = corr2(squeeze(boldtrials(i,:))',squeeze(mean(lowdiff(i,j,elecs,:),3))) ; 
    end
end
figure,imagesc(corrs) ; 

% remove outliers in both datasets ;
for i=1:3 ;  
   goodbolds{i} = find(abs(squeeze(zscore(boldtrials(i,:))))<2) ; 
end
for i=1:3 ;
    goodeegs{i} = find(abs(zscore(squeeze(mean(mean(abs(lowdiff(i,30:end,elecs,:)),2),3))))<2) ;
end
goods{1} = intersect(goodeegs{1},goodbolds{1}) ; goods{2} = intersect(goodeegs{2},goodbolds{2}) ; goods{3} = intersect(goodeegs{3},goodbolds{3}) ; 
clear corrs ; 
for i=1:3
    for j=1:size(lowdiff,2) 
        corrs(i,j) = corr2(squeeze(boldtrials(i,goods{i}))',squeeze(mean(lowdiff(i,j,elecs,goods{i}),3))) ; 
    end
end
figure,imagesc(corrs) ; 
%}