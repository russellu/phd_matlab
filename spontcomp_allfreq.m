clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie'} ; 

comps = {1:64,1:64,1:64,1:64,1:64,1:64} ; % all subjects, right and left
fcomps = {[64,80,12],[41,31,40],[61,48,15],[6,84,79],[79,80,50],[75,77,41]} ; 

latvis_fmri = [25,23,40,21,50,20] ; 
latvis_eeg = [26,14;16,6;11,23;7,22;17,23;17,7] ;

for sub=1:length(subs)

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
    prefix = 'allfreq_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ;gammas2=dir([prefix,'preproc*gamma*set']) ;
    allstims1=dir([prefix,'preproc*allstim*set']) ;allstims2=dir([prefix,'preproc*allstim*set']) ;
    movies=dir([prefix,'preproc*movie*set']) ; rests=dir([prefix,'preproc*rest*set']) ;
    setNames = {allstims1(1).name,allstims2(2).name,gammas1(1).name,gammas2(2).name,movies(1).name,rests(1).name} ;   
    
    finalersp = zeros(3,64,50,200) ; 
    clear allcorrs ; 
    for setN = 1:6;    
        EEG = pop_loadset(setNames{setN}) ; 
        if setN == 1 ; alleegs{sub} = EEG ; end

        %%% lots of noise in the data. need to remove noisy epochs in analyzer,
        %%% and then re-run ICA. especially on valerie
        trigsr = {EEG.urevent.type} ; latsr = cell2mat({EEG.urevent.latency}) ; 
        r128s = find(strcmp(trigsr,'R128')) ;
        rlats = latsr(r128s) ; 
        firstlats(setN) = rlats(1) ; 
        lastlats(setN) = rlats(length(rlats)) ; 
        
        EEG = pop_select(EEG,'nopoint',[1,rlats(1) ; rlats(length(rlats)),size(EEG.data,2)]) ; 
        EEG = eeg_checkset(EEG) ; 
    
     % get the continuous power in all frequencies
        eegsecs = EEG.pnts./EEG.srate ; 
        eegTRs = round(eegsecs./0.693) ; 
        clear epow
        for c=1:length(comps{sub})  
            [epow(c,:,:),itc,powbase,times,freqs,~,~] = newtimef(EEG.icaact(comps{sub}(c),:),EEG.pnts,[EEG.xmin,EEG.xmax],EEG.srate,0,...
                'plotersp','off','plotitc','off','baseline',NaN,'timesout',eegTRs,'nfreqs',50,'freqs',[1,100],'winsize',EEG.srate) ; 
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
           allcorrs(setN,i,:,:) = corr(squeeze(convpow(:,i,50:end-50))',squeeze(cmix(50:end-50,:))) ;        
        end
           if setN==3 || setN==4
                trigs = {'S  1','S  2','S  3'} ; 
                alltrigs = {EEG.urevent.type} ; 
                lats = cell2mat({EEG.urevent.latency}) ; 
                clear stimlats
                for t=1:length(trigs)
                    tr = find(strcmp(trigs{t},alltrigs)) ;
                    stimlats(t,:) = lats(tr) ; 
                end
                stimlats = (stimlats - firstlats(setN))./EEG.srate ; 
                stimtrs = round(stimlats./0.693) ; 
                fmricomps = segmix{setN} ; 
                fmricomps = eegfiltfft(fmricomps',1/.693,0.02,1.5) ; fmricomps = fmricomps'  ;
                % first, create the ideal time series, convolve, and correlate to
                % isolate components that were modulated by the stimulus presentation
                ideal = zeros(1,735) ; 
                stimduration = round(5./0.693) ; 
                res_stimtrs = reshape(stimtrs,[1,numel(stimtrs)]) ; 
                for i=1:length(res_stimtrs) ; ideal(res_stimtrs(i):res_stimtrs(i)+stimduration) = 1 ; end
                hrf = spm_hrf(0.693) ; 
                convedf = conv(ideal,hrf,'full') ; convedf = convedf(1:735) ;        
                corrs = corr(convedf(50:end-50)',fmricomps(50:end-50,:)) ; 
                subcorrs{sub}(setN-2,:) = corrs ; 
           end

            

        clear xc ; 
        for f=1:length(freqs)
            xc(f,:)  = crosscorr(smooth(squeeze(mean(epow(1,freqs>f & freqs<f+5,50:end-50),2)),10),squeeze(cmix(50:end-50,fcomps{sub}(1))),20) ; 
        end
        allxc(sub,setN,:,:) = xc ; 

    end

  %  subplot(2,3,sub) ; imagesc(xc) ; 
    allfinalersp(sub,:,:,:,:) = finalersp ; 
    allsubcorrs{sub} = allcorrs ; 
     
  %  figure,imagesc(squeeze(mean(epow,1)))
end

%{
fs = find(freqs>8 & freqs<25) ; 
for i=1:6 ; 
    subplot(2,3,i) ; 
    eegmat = squeeze(mean(allsubcorrs{i}(:,:,:,postcs(i)),1)) ; 
    [sv,si] = (sort(squeeze(mean(eegmat(fs,:),1)))) ; 
    topoplot(alleegs{i}.icawinv(:,si(1)),alleegs{i}.chanlocs) ; 
    latvis(i,:) = si(1:2) ; 
end
%}
for i=1:6
    plot(squeeze(mean(mean(allsubcorrs{i}(:,:,latvis_eeg(i,:),latvis_fmri(i)),1),3))) ; hold on ; 
    rels(i,:) = squeeze(mean(mean(allsubcorrs{i}(:,:,latvis_eeg(i,:),latvis_fmri(i)),1),3)) ; 
    taskrels(i,:,:) = squeeze(mean(allsubcorrs{i}(:,:,latvis_eeg(i,:),latvis_fmri(i)),3)) ; 

end

imagesc(freqs,1:6,rels,[-.35,.35]) ; xlabel('frequency(hz)') ; ylabel('subject') ; 
ylabs = {'retino#1','retino#2','gamma#1','gamma#2','movie','rest(EC)'} ; 
sublabs = {'S1','S2','S3','S4','S5','S6'} ; 
for i=1:6 ; 
    subplot(2,3,i) ; imagesc(freqs,1:6,squeeze(taskrels(i,:,:)),[-.4,.4]) ; 
    xlabel('frequency(hz)') ; set(gca,'YTick',1:6,'YTickLabel',ylabs) ; title(['subject=',num2str(i)]) ; 
end

subplot(1,2,1) ; 
imagesc(freqs,1:6,squeeze(mean(taskrels,1)),[-.3,.3]) ;  xlabel('frequency(hz)') ; set(gca,'YTick',1:6,'YTickLabel',ylabs) ; title('all tasks, mean across subjects') ; 
subplot(1,2,2) ; 
imagesc(freqs,1:6,squeeze(mean(taskrels,2)),[-.3,.3]) ; xlabel('frequency(hz)') ; set(gca,'YTick',1:6,'YTickLabel',sublabs) ; title('all subjects, mean across tasks') ; 





for i=1:length(allmeans)
    [cx,cy,cz] = centmass3(allmeans{i}) ; 
    roislice = mean(allweights{i}(:,:,cz-5:cz+5,latvis_fmri(i)),3) ; 
    subplottight(1,6,i) ; 
    plotoverlayIntensity2D(squeeze(allmeans{i}(:,:,cz-5)),mat2gray(roislice),roislice,270) ;
    title(['S',num2str(i)]) ; 
end


for i=1:length(allmeans)
    subplottight(2,6,i) ; 
    topoplot(alleegs{i}.icawinv(:,latvis_eeg(i,1)),alleegs{i}.chanlocs) ; 
    title(['S',num2str(i)]) ; 
    subplottight(2,6,i+6) ; 
    topoplot(alleegs{i}.icawinv(:,latvis_eeg(i,2)),alleegs{i}.chanlocs) ; 
end















