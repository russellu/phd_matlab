clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie'} ; 
eegcomps = {[30,26],[32,23],[21,28],[47,46],[19,37],[47,38]} ;
fcomps = {[64,80,12],[41,31,40],[61,48,15],[6,84,79],[79,80,50],[75,77,41]} ; 

tasklabs = {'retino#1 (event-related)','retino#1 (event-related)','gamma#1 (event-related)','gamma#1 (event-related)','movie (continuous viewing)','rest (eyes closed)'} ; 

for sub=1:length(subs)
    %%%% BOLD processing:
    cd(['c:/shared/badger_mri/',subs{sub},'/nii/melodic']) ; 
    mix = load('melodic_mix') ; 
    weights = load_untouch_nii('melodic_IC.nii.gz') ; 
    allweights{sub} = weights.img ; 
    meansub = load_untouch_nii('mean.nii.gz') ; 
    allmeans{sub} = meansub.img ; 
    mixinds = 1:735:735*5 ; 
    for i=1:length(mixinds) ; boldinds{i} = mixinds(i):mixinds(i)+734 ; end
    boldinds{6} = 735*5+1:735*5+1+449 ; 
    for i=1:length(boldinds) ; segmix{i} = mix(boldinds{i},:) ; end
    clear s 
    for i=1:length(segmix)
        smix = segmix{i} ; 
        [s(i,:,:),f] = spectopo(smix',0,1./0.693,'plot','off') ; 
        size(s)
    end
    figure, plot(squeeze(s(:,fcomps{sub}(1),1:100))') ; legend(tasklabs) ; 
    alls(sub,:,:) = squeeze(mean(s(:,fcomps{sub}(1:2),1:100),2))' ; 
    %%%% EEG processing:
   
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    prefix = 'highfreq_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ;gammas2=dir([prefix,'preproc*gamma*set']) ;
    allstims1=dir([prefix,'preproc*allstim*set']) ;allstims2=dir([prefix,'preproc*allstim*set']) ;
    movies=dir([prefix,'preproc*movie*set']) ; rests=dir([prefix,'preproc*rest*set']) ;

    setNames = {gammas1(1).name,gammas1(2).name,allstims1(1).name,allstims1(2).name,movies(1).name,rests(1).name} ; 
    for sn=1:length(setNames)
        EEG = pop_loadset(setNames{sn}) ; 
        [es(sub,sn,:,:),ef] = spectopo(EEG.icaact(:,30*EEG.srate:end-30*EEG.srate),0,EEG.srate,'plot','off') ;     
    end

end

for i=1:6  ;
    alles(i,:,:) = (squeeze(mean(es(i,:,eegcomps{i},:),3))') ; 
end


fs = find(ef>1 & ef<58) ; fsfreqs = ef(fs) ; 
plot(squeeze(mean(alles(:,fs,:),1)),'LineWidth',3) ;  xlim([0,length(fs)]) ; legend(tasklabs) ; set(gca,'XTick',1:10:length(fs),'XTickLabel',round(ef(1:10:end))) ; 
ylabel('log power') ; xlabel('frequency(hz)')  ; title('grand average EEG power spectrum during different tasks') ; 

%errorbar(squeeze(mean(alles(:,fs,:),1)),squeeze(std(alles(:,fs,:),0,1))./sqrt(12)) ;

plot(squeeze(mean(alls,1)),'LineWidth',3) ; legend(tasklabs) ; 
set(gca,'XTick',1:10:100,'XTickLabel',f(1:10:100)) ; title('grand average BOLD power spectrum during different tasks') ; 
ylabel('log power') ; xlabel('frequency(hz)') ; 


malles = squeeze(mean(alles,1)) ; 
malls = squeeze(mean(alls,1)) ; 


for i=1:513
    for j=1:100
        corrs(i,j) = corr2(squeeze(alles(:,i,:)),squeeze(alls(:,j,:))) ; 
    end
end

figure,
imagesc(f(1:100),ef,corrs) ; xlabel('FMRI frequency(hz)') ; ylabel('EEG frequency(hz)') ; 
imagesc(corrs) ; 




resfmri = reshape(alls(:,46,:),[1,36]) ; 
reseeg = reshape(alles(:,220,:),[1,36]) ; 
plot(resfmri,reseeg,'o') ; title(['r= ',num2str(corr2(resfmri,reseeg))]) ; 
lsline







