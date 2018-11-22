clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie'} ; 

eegcomps = {[24,25],[17],[11,15,21],[24],[27,45],[33,18]} ; 
fmricomps = {[64],[41],[61],[67],[79],[75]} ; 

for sub=1%:length(subs)

    %%%% BOLD FMRI processing:
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
    
    %%%% EEG processing:
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    prefix = '40hz_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ;gammas2=dir([prefix,'preproc*gamma*set']) ;
    allstims1=dir([prefix,'preproc*allstim*set']) ;allstims2=dir([prefix,'preproc*allstim*set']) ;
    movies=dir([prefix,'preproc*movie*set']) ; rests=dir([prefix,'preproc*rest*set']) ;
    setNames = {allstims1(1).name,allstims2(2).name,gammas1(1).name,gammas2(2).name,movies(1).name,rests(1).name} ;   
    
    finalersp = zeros(3,64,50,200) ; 
    clear allcorrs ; 
    for setN = 1:6 ;    
        EEG = pop_loadset(setNames{setN}) ; 
        if setN == 1 ; alleegs{sub} = EEG ; end
        %%% ground truth for each subject: epoch and check the gamma ERSP
        if setN==3 || setN==4
            trigs = {'S  1','S  2','S  3'} ; 
            clear ersp
            for trig=1:3
                ep = pop_epoch(EEG,{trigs{trig}},[-2,7]) ; 
                for c=1:64 
                    [ersp(trig,c,:,:),itc,powbase,times2,freqs2,~,~] = newtimef(ep.icaact(c,:),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                        'plotersp','off','plotitc','off','baseline',0,'timesout',200,'nfreqs',50,'freqs',[1,100],'winsize',64) ; 
                end
            end
            finalersp = finalersp + ersp/2 ; 
        end
        %%% lots of noise in the data. need to remove noisy epochs in analyzer,
        %%% and then re-run ICA. especially on valerie
        trigsr = {EEG.urevent.type} ; latsr = cell2mat({EEG.urevent.latency}) ; 
        r128s = find(strcmp(trigsr,'R128')) ;
        rlats = latsr(r128s) ; 
        
        firstlats(setN) = rlats(1) ; 
        
        %etrim = EEG ; %pop_select(EEG,'point',[rlats(1),rlats(length(rlats))]) ; 
        % get the continuous power in all frequencies
        eegsecs = EEG.pnts./EEG.srate ; 
        eegTRs = round(eegsecs./0.693) ; 
        clear epow
        for c=1:64 
            [epow(c,:,:),itc,powbase,times,freqs,~,~] = newtimef(EEG.icaact(c,:),EEG.pnts,[EEG.xmin,EEG.xmax],EEG.srate,0,...
                'plotersp','off','plotitc','off','baseline',NaN,'timesout',eegTRs,'nfreqs',100,'freqs',[1,100],'winsize',EEG.srate) ; 
        end   
        cmix = segmix{setN}(1:eegTRs,:) ;  
        cmix = eegfiltfft(cmix',1/.693,.02,1.5) ; cmix = cmix' ; 
        hrf = spm_hrf(0.693) ; 
        clear convpow
        for i=1:size(epow,1)
            for j=1:size(epow,2)
                conved = conv(squeeze(epow(i,j,:)),hrf,'full') ;
                convpow(i,j,:) = conved(1:eegTRs) ;             
            end
        end  
        for i=1:size(convpow,2)
           allcorrs(setN,i,:,:) = corr(squeeze(convpow(:,i,25:end-25))',squeeze(cmix(25:end-25,:))) ;        
        end
    end
    allfinalersp(sub,:,:,:,:) = finalersp ; 
    allsubcorrs{sub} = allcorrs ; 
    %figure,for j=1:64 ; subplot(5,13,j) ; imagesc(squeeze(mean(finalersp(1:2,j,:,:),1)),[-6,6]) ; title(j) ; end
end


%{
for i=1:6 ; subplot(2,6,i) ; 
    plot(squeeze(mean(allfinalersp(i,:,eegcomps{i}(1),:,50:150),5))','LineWidth',2) ; 
    xlim([0,50]) ; ylim([-8,8]) ; hline(0,'k') ; 
    set(gca,'XTick',1:5:50,'XTickLabel',round(freqs2(1:5:50))) ; xlabel('frequency(hz)') ; ylabel('power(db)') ; 
    if i==1 ; legend({'0%rnd','10%rnd','100%rnd'}) ; end   
    title(['subject',num2str(i)]) 
end
for i=1:6 ; subplot(2,6,i+6) ; 
    imagesc(times2,freqs2,squeeze(mean(allfinalersp(i,:,eegcomps{i}(1),:,:),2)),[-5,5]) ; axis xy
    xlabel('time(s)') ; ylabel('frequency(hz)') ; vline([0,5],'k') ; 
end
for i=1:6 ; 
   subplot(1,6,i) ;
   imagesc(squeeze(allsubcorrs{i}(:,:,eegcomps{i}(1),fmricomps{i}(1))),[-.5,.5])
   fcorrs(i,:,:) = squeeze(allsubcorrs{i}(:,:,eegcomps{i}(1),fmricomps{i}(1))) ; 
   hline(0.5:1:6,'k') ; 
end
%plot(squeeze(mean(mean(fcorrs(:,1:2,:)))),'r','LineWidth',2) ; hold on ; 
plot(squeeze(mean(mean(fcorrs(:,3:4,:)))),'m','LineWidth',2) ; hold on ; 
plot(squeeze(mean(fcorrs(:,5,:))),'b','LineWidth',2) ;
plot(squeeze(mean(fcorrs(:,6,:))),'k','LineWidth',2) ;
legend({'event-related','movie','rest'}) ; 
set(gca,'XTick',1:10:100,'XTickLabel',freqs(1:10:100)) ; xlabel('frequency(hz)') ; ylabel('mean correlation (r)') ; hline(0,'k') ; 
for i=1:length(allmeans)
    [cx,cy,cz] = centmass3(allmeans{i}) ; 
    roislice = mean(allweights{i}(:,:,cz-6:cz-3,fmricomps{i}(1)),3) ; 
    subplot(1,6,i) ; 
    plotoverlayIntensity2D(squeeze(allmeans{i}(:,:,cz-5)),mat2gray(roislice),roislice,270) ;

end
%}






