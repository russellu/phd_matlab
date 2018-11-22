clear all ; close all ; 
fmrivols = {'bp_mc_retino_gamma_01.nii','bp_mc_retino_gamma_02.nii'} ; 
stimtrigs = {'S  1','S  2','S  3'} ; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie'} ;
subcomps = {[28,24,20,32],[21,47],[22,23,52],[25,44],[27,26,29],[19,26,21],[44,21]} ;
name = 'dina' ;

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
%f1 = load_untouch_nii('f_mc_Test_Russell_2015_11_12_WIP_EEG-fMRI_MB3_3.75mm_SENSE_6_1.nii.gz') ;
%for i=1:size(tres,1) ; f1.img = squeeze(tres(i,:,:,:)) ; save_untouch_nii(f1,['tvals_',num2str(i),'.nii.gz']) ; end
%for i=1:size(pres,1) ; f1.img = squeeze(pres(i,:,:,:)<.0001) ; save_untouch_nii(f1,['pvals_',num2str(i),'.nii.gz']) ; end
%{
clear ersp ; 
for eset=1:length(eegsets) 
    cd(['c:/shared/badger_eeg/',name,'/']) ; 
    eeg = pop_loadset(eegsets{eset}) ; 
    for st=1:3
        ep = pop_epoch(eeg,{stimtrigs{st}},[-2,7]) ; 
        for comp=1:size(ep.icaact,1) ; 
            for trial=1:size(ep.icaact,3) 
                [ersp(st,comp,(eset*size(ep.data,3)-size(ep.data,3)) + trial,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(comp,:,trial)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                    'plotitc','off','plotersp','off','baseline',NaN,'freqs',[1,128],'nfreqs',64,'winsize',64 );         
            end
        end
    end
end
bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,200]) ; 

comps = [20,24] ; allc = zeros(1,64) ; allc(comps) = 1 ; 
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(times,freqs,squeeze(mean(mean(bersp(1:2,i,:,:,:),3),1)),[-5,5]) ; axis xy ; title(i) ; end
%lnames = {'0%rnd','10%rnd','100%rnd','mean'} ;
%for i=1:3 ; subplot(1,3,i) ; imagesc(times,freqs,squeeze(mean(mean(bersp(i,comps,:,:,:),2),3)),[-6,6]) ; axis xy ; title(lnames{i}) ; xlabel('time(s)' ) ;ylabel('frequency(hz)') ;  end ; suptitle('mean ERSP') ; 

%%% get the subtracted power values
elecs = [15,51,7,37,19,38,8,52,16,59,45,31,46,60,9,20,10] ; es = zeros(64,1) ; es(elecs) = 1 ; 
clear baselow tasklow allepochs
for eset=1:length(eegsets) 
    cd(['c:/shared/badger_eeg/',name,'/']) ; 
    eeg = pop_loadset(eegsets{eset}) ; 
    eeg = pop_subcomp(eeg,find(allc==0)) ; 
    fcount = 1 ; 
    for f=1:120
    eeglow = eeg ; eeglow.data = eegfiltfft(eeglow.data,eeglow.srate,f,f+6) ; 
    for st=1:3
        ep = pop_epoch(eeglow,{stimtrigs{st}},[-6,5]) ; 
        lowdata = ep.data.^2 ; 
        baselow(st,fcount,:,(eset*16-16)+1:eset*16) = squeeze(mean(lowdata(:,ep.times<0,:),2)) ; 
        tasklow(st,fcount,:,(eset*16-16)+1:eset*16) = squeeze(mean(lowdata(:,(ep.times>0) & ep.times<5000,:),2)) ; 
        %allepochs(st,fcount,:,:,(eset*16-16)+1:eset*16) = squeeze(lowdata(elecs,:,:)) ; 
    end
    fcount = fcount + 1 ; 
    end
end
lowdiff = tasklow-baselow ; 
%mepochs= squeeze(mean(allepochs,3)) ; 

%plot(squeeze(mean(mepochs(1,30:60,:,:),2))) ; 
%plot(mean(imfilter(squeeze(mean(mepochs(1,30:50,:,:),2)),fspecial('gaussian',[9,1],3)),2)) ;
%}

comps = subcomps{3} ; allc = zeros(1,64) ; allc(comps) = 1 ;
elecs = [15,51,7,37,19,38,8,52,16,59,45,31,46,60,9,20,10] ; es = zeros(64,1) ; es(elecs) = 1 ; 
for eset=1:length(eegsets) 
    cd(['c:/shared/badger_eeg/',name,'/']) ; 
    eeg = pop_loadset(eegsets{eset}) ; 
    eeg = pop_subcomp(eeg,find(allc==0)) ;      
    for st=1:3
        ep = pop_epoch(eeg,{stimtrigs{st}},[-2,7]) ;
        for elec=1:length(elecs)
            for trial=1:size(ep.data,3)
                [ersp(st,elec,(eset*size(ep.data,3)-size(ep.data,3)) + trial,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.data(elecs(elec),:,trial)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                    'plotitc','off','plotersp','off','baseline',NaN,'freqs',[1,120],'nfreqs',64,'winsize',64 ) ;
            end
        end
    end
end
bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,200]) ; 


tthresh=7 ; 
roi = ((squeeze(mean(tres,1))>tthresh)) ; roivox = find(roi==1) ; 
for i=1:size(taskvols,1)
    for j=1:size(taskvols,5)
        boldij = squeeze(diffvols(i,:,:,:,j)) ; 
        boldtrials(i,j) = squeeze(mean(boldij(roivox))) ; 
    end
end

mersp = squeeze(mean(bersp,2)) ; 
mtersp = squeeze(mean(mersp(:,:,:,times>0 & times<5),4)) ; 
clear corrs
for i=1:3
    for j=1:64 ;
        for k=1:200
            corrs(i,j,k) = corr2(squeeze(mersp(i,:,j,k)),boldtrials(i,:)) ; 
        end
    end
end

% look at all trials at once
clear allbold allmt ; 
allbold(1:32) = boldtrials(1,:) ; allbold(33:64) = boldtrials(2,:) ; allbold(65:96) = boldtrials(3,:) ; 
allmt(1:32,:) = squeeze(mtersp(1,:,:)) ; allmt(33:64,:) = squeeze(mtersp(2,:,:)) ; allmt(65:96,:) = squeeze(mtersp(3,:,:)) ; 
[sv,si] = sort(allbold,'descend') ; 
for i=1:64 ; subplot(5,13,i) ; bar(smooth(allmt(si,i))) ; end


for i=1:3
   [~,si] = sort(boldtrials(i,:),'descend') ;  
    sortersp(i,:,:) = mtersp(i,si,:) ; 
end
bar(squeeze(mean(abs(sortersp(:,:,:)),3)))


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

%{
mbersp  = squeeze(mean(bersp(:,:,:,:,times>0 & times<5),5)) ; 
for i=1:3 ; 
    for j=1:64 ; 
        for k=1:64  
            corrs(i,j,k) = corr2(squeeze(mbersp(i,j,:,k)),boldtrials(i,:)') ; 
        end
    end
end
%}

clear corrs ; 
for i=1:3
    for j=1:size(lowdiff,2) 
        corrs(i,j) = corr2(squeeze(boldtrials(i,:))',squeeze(mean(lowdiff(i,j,elecs,:),3))) ; 
    end
end
figure,imagesc(corrs) ; 

figure,subplot(1,2,1) ; bar(boldtrials)
subplot(1,2,2) ; bar(squeeze(mean(mean(lowdiff(:,40:60,elecs,:),2),3)))

% correlate with ERSP:
for c=1:64 ; 
mbersp  = squeeze(mean(bersp(:,c,:,:,:),2)) ; 
for i=1:3 ; 
    for j=1:64 ; 
        for k=1:200
            bcorrs(i,j,k) = corr2(squeeze(mbersp(i,:,j,k)),squeeze(boldtrials(i,:))) ; 
        end
    end
end
figure,
for i=1:3 ; subplot(2,2,i) ; imagesc(squeeze(bcorrs(i,:,:)),[-1,1]) ; end
suptitle(num2str(c)) ; 
end


figure, 
subplot(1,2,1) ;
hold on ; plot(smooth(corrs(1,:)),'b') ; plot(smooth(corrs(2,:)),'g') ; plot(smooth(corrs(3,:)),'r') ; hline(0,'k') ; 
subplot(1,2,2) ; plot(squeeze(mean(mean(mean(bersp(:,comps,:,1:60,times>0 & times<5),2),3),5))') ; hline(0,'k') ; 

meaneeg = squeeze(mean(lowdiff(:,:,elecs,:),3)) ; 
meancorrs = squeeze(mean(corrs,1)) ; 
for i=1:3 
   hrp(i,:) = sum(squeeze(meaneeg(i,:,:)).*repmat(meancorrs,[32,1])') ;  
end


% sort and visualize
for i=1:3 ; [~,si(i,:)] = sort(boldtrials(i,:),'descend') ; end ; 
for i=1:5:120 ; 
    figure ; 
    subplot(2,2,1) ; bar(squeeze(mean(meaneeg(1,i:i+5,si(1,:)),2))) ; 
    subplot(2,2,2) ; bar(squeeze(mean(meaneeg(2,i:i+5,si(2,:)),2))) ; 
    subplot(2,2,3) ; bar(squeeze(mean(meaneeg(3,i:i+5,si(3,:)),2))) ; 
end


plot(corrs','LineWidth',2) ; hline(0,'k') ; hold on ; plot(mean(corrs),'LineWidth',3,'Color',[0,0,0]) ; legend(lnames) ; xlabel('frequency(hz)') ; ylabel('correlation') ; 
title('EEG-BOLD relationship, single trial correlations in 1 subject') ; 




for i=1:5:120 ; 
figure ; 
plot(squeeze(mean(meaneeg(1,i:i+5,:),2)),boldtrials(1,:),'.') ; hold on ; 
plot(squeeze(mean(meaneeg(2,i:i+5,:),2)),boldtrials(2,:),'.','Color',[0,1,0]) ;
plot(squeeze(mean(meaneeg(3,i:i+5,:),2)),boldtrials(3,:),'.','Color',[1,0,0]) ;
title(i)  ;
end

topoplot(zeros(1,64),eeg.chanlocs,'electrodes','numbers') ; 
elecs = [59,45,31,46,60,9,20,10] ; 
meantask = squeeze(mean(taskhigh(:,elecs,:),2)) ; 
plot(meantask(2,:),boldtrials(2,:),'o') ; 
for i=1:65 ; 
    subplot(7,10,i) ; 
    reseeg = reshape(mean(mean(lowdiff(:,i:i+5,elecs,:),2),3),[1,32*3]) ; 
    resfmri = reshape(boldtrials,[1,32*3]) ; 
    plot(resfmri,reseeg,'.') ; title(num2str(corr2(resfmri,reseeg))) ; 
end


subplot(1,2,1) ; plot(squeeze(mean(mean(lowdiff(:,5:13,elecs,:),2),3))',boldtrials(:,:)','o','LineWidth',2) ; legend(lnames(1:3)) ; 
eegc = reshape(squeeze(mean(mean(lowdiff(:,5:13,elecs,:),2),3)),[1,32*3]) ; boldc = reshape(boldtrials,[1,32*3]) ; 
title(['low frequencies (5-13Hz), corr2 = ',num2str(corr2(eegc,boldc))]) ; 

xlabel('EEG summed power (task-rest)') ; ylabel('BOLD task-rest') ; 
subplot(1,2,2) ; plot(squeeze(mean(mean(lowdiff(:,30:60,elecs,:),2),3))',boldtrials(:,:)','o','LineWidth',2) ; legend(lnames(1:3)) ; 
eegc = reshape(squeeze(mean(mean(lowdiff(:,30:60,elecs,:),2),3)),[1,32*3]) ; boldc = reshape(boldtrials,[1,32*3]) ; 
title(['high frequencies (30-60Hz), corr2 = ',num2str(corr2(eegc,boldc))]) ; 

