clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 

comps = {[30,26],[32,23],[21,28],[47,46],[19,37],[47,38],[44,46]} ; % all subjects, right and left
fcomps = {[64,80,12],[41,31,40],[61,48,15],[6,84,79],[79,80,50],[75,77,41],[48,50]} ; 

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
    cd .. 
    movie = load_untouch_nii('bp_*_retino_movie.nii.gz') ; 
    rest = load_untouch_nii('bp_*_retino_rest.nii.gz') ; 

    %%%% EEG processing:
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    prefix = 'highfreq_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ;gammas2=dir([prefix,'preproc*gamma*set']) ;
    allstims1=dir([prefix,'preproc*allstim*set']) ;allstims2=dir([prefix,'preproc*allstim*set']) ;
    movies=dir([prefix,'preproc*movie*set']) ; rests=dir([prefix,'preproc*rest*set']) ;
    setNames = {allstims1(1).name,allstims1(2).name,gammas1(1).name,gammas1(2).name,movies(1).name,rests(1).name} ;   

    hrf = spm_hrf(0.693) ; 
    
    for setN = 1:6;    
        EEG = pop_loadset(setNames{setN}) ; 
        lats = {EEG.urevent.latency} ; 
        types = {EEG.urevent.type} ; 
        r128s = find(strcmp('R128',types)) ; 
        rlats = cell2mat(lats(r128s)) ; 
        
        eegacts = EEG.icaact(:,rlats(1):rlats(length(rlats))) ; 
        eegsecs = size(eegacts,2)./EEG.srate ; 
        eegTRs = round(eegsecs./0.693) ; 
    
        clear epow
        for c=1:length(comps{sub})  
            [epow(c,:,:),itc,powbase,times,freqs,~,~] = newtimef(eegacts(comps{sub}(c),:),size(eegacts,2),[0,size(eegacts,2)./EEG.srate],EEG.srate,0,...
                'plotersp','off','plotitc','off','baseline',NaN,'timesout',eegTRs,'nfreqs',50,'winsize',250,'freqs',[1,100]) ; 
        end
        indlags = -20:20 ; 
        mepow = squeeze(mean(epow,1)) ; 
        smix = segmix{setN} ; smix = eegfiltfft(smix',1/0.693,0.03,1.5) ; smix = smix' ; 
        smix = mean(smix(1:size(mepow,2),fcomps{sub}(1:2)),2) ; 
        clear restcorrs
        hrf = spm_hrf(0.693) ; 
        if setN == 6 
            restmix(sub,:) = smix ; 
            restpow(sub,:,:) = mepow ; 
            restimg = rest.img(:,:,:,1:size(mepow,2)) ; 
            for i=1:50 ; 
                powi = (squeeze(mean(mepow(freqs>freqs(i)-2.5 & freqs<freqs(i)+2.5,:),1))) ; 
                powi2 = conv(powi,hrf,'full') ; powi = powi2(1:length(powi)) ; 
                restcorrs(:,:,:,i) = voxcorr(restimg(:,:,:,40:end-40),powi(40:end-40)) ; 
            end
            
        end 
        if setN == 5  
            moviemix(sub,:) = smix ; 
            moviepow(sub,:,:) = mepow ; 
            movieimg = movie.img(:,:,:,1:size(mepow,2)) ; 
            for i=1:50 ; 
                powi = (squeeze(mean(mepow(freqs>freqs(i)-2.5 & freqs<freqs(i)+2.5,:),1))) ; 
                powi2 = conv(powi,hrf,'full') ; powi = powi2(1:length(powi)) ; 
                moviecorrs(:,:,:,i) = voxcorr(movieimg(:,:,:,40:end-40),powi(40:end-40)) ; 
            end
            
            
        end 
        for i=1:length(indlags) % for all time lags
            for j=1:size(mepow,1) % for all freqs
                powj = smooth(mean(mepow(freqs>freqs(j)-2.5 & freqs<freqs(j)+2.5,:))) ; 
                corrs(j,i) = corr2(powj(50:end-50)',smix(50+indlags(i):end-50+indlags(i))') ;               
            end   
        end
        allcorrs(sub,setN,:,:) = corrs ; 
    end
    figure,for i=1:50  ; subplot(5,10,i) ; imagesc(squeeze(mean(restcorrs(:,:,8:15,i),3)),[-.25,.25]) ; end
    allmovies{sub} = moviecorrs ; 
    allrests{sub} = restcorrs ; 
end

for i=1:7 ; figure
    for j=1:size(allrests{i},3)
        subplot(4,10,j) ; imagesc(squeeze(mean(allrests{i}(:,:,j,15:25),4)),[-.2,.2]) ; 
        
    end
end


for sub=1:length(subs) ; 
cd(['c:/shared/badger_mri/',subs{sub},'/nii/']) ; 
freqref = load_untouch_nii('freq_ref.nii.gz') ; 
freqref.img = allmovies{sub} ; 
save_untouch_nii(freqref,'moviefreqcorrs.nii.gz') ; 

end

















