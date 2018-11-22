clear all ; close all ; 
cc = [19,20,60] ; 
eegsets = {'alex_retino_gamma_01.set','retino_gamma_02.set'} ; 
fmrivols = {'bp_reg_mc_Test_Russell_2015_11_12_WIP_EEG-fMRI_MB3_3.75mm_SENSE_6_1.nii.gz','bp_reg_mc_Test_Russell_2015_11_12_WIP_EEG-fMRI_MB3_3.75mm_SENSE_8_1.nii.gz'} ; 
cd ('c:/shared/badger/alex/Russell BADGER 2015-11-12/mel_gamma') ; mix = load('melodic_mix') ; mixes(1,:) = mix(1:735,4) ; mixes(2,:) = mix(736:end,4) ; 
clear corrvols diffvols ; 
stimtrigs = {'S  1','S  2','S  3'} ; 
stimtrigcounts = [1,1,1] ; 
for eset = 1:length(eegsets) 
cd c:/shared/badger/alex/ ; ls 
eeg = pop_loadset(eegsets{eset}) ; 
cd('Russell BADGER 2015-11-12') ; 
bold = load_untouch_nii(fmrivols{eset}) ; 
eegevents = {eeg.urevent.type} ; 
eeglats = {eeg.urevent.latency} ; 
volonsets = find(strcmp('R128',eegevents)) ; 
volonsetlats = cell2mat(eeglats(volonsets)) ; 
basetime = 2 ; baseTR = ceil(basetime/bold.hdr.dime.pixdim(5)) ; 
tasktime = 12 ; taskTR = ceil(tasktime/bold.hdr.dime.pixdim(5)) ; 
hrf = spm_hrf(bold.hdr.dime.pixdim(5)) ; 
task = zeros(1,taskTR) ; task(1:round(5/bold.hdr.dime.pixdim(5))) = 1 ; task = conv(task,hrf) ; task = task(1:taskTR) ; 
peakinds = find(task==max(task))-3:find(task==max(task))+3 ; 

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
       mixvals(st,:,stimtrigcounts(st)) = squeeze(mixes(eset,mininds(i):mininds(i)+taskTR-1)) ;
       stimtrigcounts(st) = stimtrigcounts(st) + 1 ; 
    end
end

end

restask = reshape(taskvols,[numel(taskvols(:,:,:,:,1)),32]) ; 
resbase = reshape(basevols,[numel(basevols(:,:,:,:,1)),32]) ; 
[h,p,ci,stats] = ttest(restask',resbase') ; 

f1 = load_untouch_nii('f_mc_Test_Russell_2015_11_12_WIP_EEG-fMRI_MB3_3.75mm_SENSE_6_1.nii.gz') ;
tstats = stats.tstat ; 
tres = reshape(tstats,[3,size(corrvols,2),size(corrvols,3),size(corrvols,4)]) ; 
pres = reshape(p,[3,size(corrvols,2),size(corrvols,3),size(corrvols,4)]) ;
for i=1:size(tres,1) ; f1.img = squeeze(tres(i,:,:,:)) ; save_untouch_nii(f1,['tvals_',num2str(i),'.nii.gz']) ; end
for i=1:size(pres,1) ; f1.img = squeeze(pres(i,:,:,:)<.0001) ; save_untouch_nii(f1,['pvals_',num2str(i),'.nii.gz']) ; end


%%%% 
clear ersp ; 
for eset=1:length(eegsets) 
    cd c:/shared/badger/alex/ ; 
    eeg = pop_loadset(eegsets{eset}) ; 
    for st=1:3
        ep = pop_epoch(eeg,{stimtrigs{st}},[-7,7]) ; 
        for comp=1:size(ep.icaact,1) ; 
            for trial=1:size(ep.icaact,3) 
                [ersp(st,comp,(eset*size(ep.data,3)-size(ep.data,3)) + trial,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(comp,:,trial)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                    'plotitc','off','plotersp','off','baseline',NaN,'freqs',[1,128],'nfreqs',64,'winsize',64 );         
            end
        end
    end
end
bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,200]) ; 

comps = [18,27,28] ; allc = zeros(1,64) ; allc(comps) = 1 ; 
lnames = {'0%rnd','10%rnd','100%rnd','mean'} ;
for i=1:3 ; subplot(1,3,i) ; imagesc(times,freqs,squeeze(mean(mean(bersp(i,comps,:,:,:),2),3)),[-6,6]) ; axis xy ; title(lnames{i}) ; xlabel('time(s)' ) ;ylabel('frequency(hz)') ;  end ; suptitle('mean ERSP') ; 

%%% get the subtracted power values
clear baselow tasklow
for eset=1:length(eegsets) 
    cd c:/shared/badger/alex/ ; 
    eeg = pop_loadset(eegsets{eset}) ; 
    eeg = pop_subcomp(eeg,find(allc==0)) ; 
    fcount = 1 ; 
    for f=6:85
    eeglow = eeg ; eeglow.data = eegfiltfft(eeglow.data,eeglow.srate,f-5,f+5) ; 
    for st=1:3
        ep = pop_epoch(eeglow,{stimtrigs{st}},[-5,6]) ; 
        lowdata = ep.data.^2 ; 
        baselow(st,fcount,:,(eset*16-16)+1:eset*16) = squeeze(mean(lowdata(:,find(ep.times<0),:),2)) ; 
        tasklow(st,fcount,:,(eset*16-16)+1:eset*16) = squeeze(mean(lowdata(:,find(ep.times>0),:),2)) ; 
    end
    fcount = fcount + 1 ; 
    end
end
lowdiff = tasklow-baselow ; 

tthresh=8 ; 
roi = ((squeeze(mean(tres,1))>tthresh)) ; roivox = find(roi==1) ; 
for i=1:size(taskvols,1)
    for j=1:size(taskvols,5)
        boldij = squeeze(diffvols(i,:,:,:,j)) ; 
        boldtrials(i,j) = squeeze(mean(boldij(roivox))) ; 
    end
end

clear corrs ; 
for i=1:3
    for j=1:size(lowdiff,2) 
        corrs(i,j) = corr2(squeeze(boldtrials(i,:))',squeeze(mean(lowdiff(i,j,elecs,:),3))) ; 
    end
end
imagesc(corrs)
plot(corrs','LineWidth',2) ; hline(0,'k') ; hold on ; plot(mean(corrs),'LineWidth',3,'Color',[0,0,0]) ; legend(lnames) ; xlabel('frequency(hz)') ; ylabel('correlation') ; 
title('EEG-BOLD relationship, single trial correlations in 1 subject') ; 

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



%{
rndpercs = [0,10,100] ; rndnames = {'0%rnd','10%rnd','100%rnd'} ; 
mbersp = squeeze(mean(bersp(:,comps,:,:,:),2)) ; 
for i=1:3 ; 
   subplot(1,3,i) ; 
   plot(squeeze(mean(mean(mbersp(i,:,freqs>10 & freqs<25,times>0 & times<5),3),4)), squeeze(mean(mean(mbersp(i,:,freqs>35 & freqs<70,times>0 & times<5),3),4)),'o') ; 
   title(['corr2 = ',num2str(corr2(squeeze(mean(mean(mbersp(i,:,freqs>10 & freqs<25,times>0 & times<5),3),4)), squeeze(mean(mean(mbersp(i,:,freqs>35 & freqs<70,times>0 & times<5),3),4))))]) ;
   xlabel('low(10-25Hz) frequencies') ; ylabel('high(35-70Hz) frequencies') ; 
end

errorbar(squeeze(mean(mean(mbersp(:,:,:,times>0 & times<5),2),4))',squeeze(std(mean(mbersp(:,:,:,times>0 & times<5),4),0,2))'./sqrt(32),'LineWidth',3) ; legend(rndnames) ; hline(0,'k') ; 
xlim([1,64]) ; xlabel('frequency(hz)') ; ylabel('power(db)') ; set(gca,'XTick',1:4:64,'XTickLabel',round(freqs(1:4:64))) ;

%%%%%% compare EEG and BOLD
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(bersp(1,i,:,:,:),3)),[-6,6]) ; title(i) ; end

for comp=1:64
mbersp = squeeze(mean(bersp(:,comp,:,:,:),2)) ; 
tthresh = 5 ; 
roivox = find((squeeze(mean(tres,1))>tthresh)) ; roi = ((squeeze(mean(tres,1))>tthresh)) ; roivox = find(roi==1) ; 
for i=1:size(diffvols,1)   
    for j=1:size(diffvols,5)
        volij = squeeze(diffvols(i,:,:,:,j)) ; 
        roidiff(i,j) = squeeze(mean(volij(roivox))) ;      
    end
end
meanmix = squeeze(mean(mixvals(:,9:12,:),2)) ; 
for i=1:3 
    for j=1:64 ; 
        for k=1:200
            kernel(i,j,k) = corr2(squeeze(mbersp(i,:,j,k)),squeeze(meanmix(i,:))) ;
        end
    end
end
figure,for i=1:3 ; subplot(1,3,i) ; imagesc(squeeze(kernel(i,:,:))) ; end ; suptitle(num2str(comp)) ; 
%figure,imagesc(squeeze(mean(kernel(1:3,:,:))))
end

meanmix = squeeze(mean(mixvals(:,9:12,:),2))-squeeze(mean(mixvals(:,3:4,:),2)) ; 
for i=1:3 
    for j=1:64 ; 
        for k=1:200
            kernel(i,j,k) = corr2(squeeze(mbersp(i,:,j,k)),squeeze(meanmix(i,:))) ;
        end
    end
end
for i=1:3 ; subplot(1,3,i) ; imagesc(squeeze(kernel(i,:,:))) ; end
figure,imagesc(squeeze(mean(kernel(1:3,:,:))))

mtbersp = squeeze(mean(mbersp(:,:,:,times>0 & times<5),4)) ; 
for i=1:3
    for j=1:64 ;
        corrs(i,j) = corr2(squeeze(mtbersp(i,:,j)),roidiff(i,:)) ;       
    end
end


cs = [7,18,27,28] ; 
elabs = {eeg.chanlocs.labels} ; 

clear stask sbase
% get the asynchronous broadband component
for eset=1:length(eegsets) 
    cd c:/shared/badger/alex/ ; 
    eeg = pop_loadset(eegsets{eset}) ; 
    allc = zeros(1,64) ; allc(cs) = 1 ; 
    eeg = pop_subcomp(eeg,find(allc==0)) ; 
    for st=1:3
        ep = pop_epoch(eeg,{stimtrigs{st}},[-5,5]) ; 
        for comp=1:64 ; 
            [sbase(st,comp,(eset*16-16)+1:eset*16,:),fhz] = spectopo(squeeze(ep.data(comp,1:ep.pnts/2,:))',0,ep.srate,'plot','off') ;  
            [stask(st,comp,(eset*16-16)+1:eset*16,:),~] = spectopo(squeeze(ep.data(comp,ep.pnts/2:end,:))',0,ep.srate,'plot','off') ;                   
        end
    end
end

sdiff = stask-sbase ; 
imagesc(squeeze(sdiff(1,28,:,:)),[-15,15])
topoplot(squeeze(mean(mean(sdiff(2,:,:,find(fhz>30 & fhz<50)),3),4)),eeg.chanlocs) ;

%meanmix = squeeze(mean(mixvals(:,9:12,:),2))-squeeze(mean(mixvals(:,2:3,:),2)) ; 
clear corrs ; 
for i=1:64
    for f=1:256-15
    %mbersp = squeeze(mean(bersp(:,i,:,:,:),2)) ; 
    %broadband = squeeze((mean(mean(abs(mbersp(:,:,:,times>0 & times<5)),3),4))) ; 
    %broadband = squeeze(vars(:,i,:)) ; 
    broadband = squeeze(mean(sdiff(:,i,:,f:f+10),4)) ; 
    corrs(i,f,1) = corr2(roidiff(1,:),broadband(1,:)) ; 
    corrs(i,f,2) = corr2(roidiff(2,:),broadband(2,:)) ; 
    corrs(i,f,3) = corr2(roidiff(3,:),broadband(3,:)) ; 
    end
end
%cs = [7,18,27,28] ; 
figure,plot(squeeze(mean(mean(corrs(cs,:,1:2),1),3)),'LineWidth',4) ; hold on ; plot(squeeze(mean(mean(corrs(cs,:,3),1),3)),'LineWidth',4,'Color',[1,0,0]) 
set(gca,'XTick',1:8:length(fhz),'XTickLabel',round(fhz(1:8:length(fhz)))) ; hline(0,'k') ; 

%}
