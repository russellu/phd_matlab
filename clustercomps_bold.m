clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 

categories = {'bcg','center','eyes','rim'... 
              'motor_left','motor_center','motor_right',... 
              'frontal_left','frontal_center','frontal_right',... 
              'occ_left','occ_center','occ_right',... 
              'parietal_left','parietal_center','parietal_right',... 
              'temporal_left','temporal_right'... 
              'noise_symmetric','noise_focal','noise_broad'} ; 
       
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
    %movie = load_untouch_nii('warped_movie.nii.gz') ; 
    %rest = load_untouch_nii('bp_*_retino_rest.nii.gz') ; 
    rest = load_untouch_nii('antsf/warped_resting.nii.gz') ; 

    %%%% EEG processing:
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    prefix = 'allfreq_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ;gammas2=dir([prefix,'preproc*gamma*set']) ;
    allstims1=dir([prefix,'preproc*allstim*set']) ;allstims2=dir([prefix,'preproc*allstim*set']) ;
    movies=dir([prefix,'preproc*movie*set']) ; rests=dir([prefix,'preproc*rest*set']) ;
    setNames = {allstims1(1).name,allstims1(2).name,gammas1(1).name,gammas1(2).name,movies(1).name,rests(1).name} ;   

    hrf = spm_hrf(0.693) ; 
    
    for setN = 6    
        
        compinds = load('compinds') ; compinds = compinds.inds_s ; 
        
        EEG = pop_loadset(setNames{setN}) ; 
        lats = {EEG.urevent.latency} ; 
        types = {EEG.urevent.type} ; 
        r128s = find(strcmp('R128',types)) ; 
        rlats = cell2mat(lats(r128s)) ; 
        
        eegacts = EEG.icaact(:,rlats(1):rlats(length(rlats))) ; 
        eegsecs = size(eegacts,2)./EEG.srate ; 
        eegTRs = round(eegsecs./0.693) ; 
        
        clear epow
        for c=1:length(compinds)  
            [epow(c,:,:),itc,powbase,times,freqs,~,~] = newtimef(eegacts(c,:),size(eegacts,2),[0,size(eegacts,2)./EEG.srate],EEG.srate,0,...
                'plotersp','off','plotitc','off','baseline',NaN,'timesout',eegTRs,'nfreqs',50,'winsize',250,'freqs',[1,100]) ; 
        end
        indlags = -20:20 ; 
        occs = find(compinds==20| compinds==20| compinds==20) ; 
        mepow = squeeze(mean(epow(occs,:,:),1)) ; 
        %smix = segmix{setN} ; smix = eegfiltfft(smix',1/0.693,0.03,1.5) ; smix = smix' ; 
        %smix = mean(smix(1:size(mepow,2),fcomps{sub}(1:2)),2) ; 
        %{
        clear restcorrs
        hrf = spm_hrf(0.693) ; 
        if setN == 6 
           % restmix(sub,:) = smix ; 
            restpow(sub,:,:) = mepow ; 
            restimg = rest.img(:,:,:,1:size(mepow,2)) ; 
            for i=1:50 ; 
                powi = (squeeze(mean(mepow(freqs>freqs(i)-2.5 & freqs<freqs(i)+2.5,:),1))) ; 
                powi2 = conv(powi,hrf,'full') ; powi = powi2(1:length(powi)) ; 
                restcorrs(:,:,:,i) = voxcorr(restimg(:,:,:,40:end-40),powi(40:end-40)) ; 
            end
        end 
        %}
        % do the alpha band power for all components 
        clear restcorrs
        occs = find(compinds>10 & compinds<13) ; 
        mepow = squeeze(mean(epow(occs,:,:),1)) ; 
        hrf = spm_hrf(0.693) ; 
        restimg = rest.img(:,:,:,1:size(mepow,2)) ; 
        for i=1:50 ; 
            powi = (squeeze(mean(mepow(freqs>freqs(i)-2.5 & freqs<freqs(i)+2.5,:),1))) ; 
            powi2 = conv(powi,hrf,'full') ; powi = powi2(1:length(powi)) ; 
            restcorrs(:,:,:,i) = voxcorr(restimg(:,:,:,40:end-40),powi(40:end-40)) ; 
        end

        %{
        if setN == 5  
           % moviemix(sub,:) = smix ; 
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
          %}
    end    
   % figure,for i=1:50  ; subplot(5,10,i) ; imagesc(squeeze(mean(restcorrs(:,:,8:15,i),3)),[-.25,.25]) ; end
   % allmovies{sub} = moviecorrs ; 
    allrests{sub} = restcorrs ; 
end

%figure,for i=1:64 ; subplot(5,13,i) ; topoplot(squeeze(EEG.icawinv(:,i)),EEG.chanlocs) ; title(i) ; end
%figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(restcorrs(:,:,7:11,i),3)),[-.25,.25]) ; title(i) ; end


avgrests = zeros(7,64,64,33,50) ; 
for i=1:length(allrests)
    avgrests(i,:,:,:,:) = allrests{i} ; 
  %  avgmovies(i,:,:,:,:) = allmovies{i} ; 
end
ref = load_untouch_nii('c:/shared/badger_mri/jeremie/nii/freqref.nii.gz') ; 
mrests = squeeze(mean(avgrests,1)) ; mrests(isnan(mrests)) = 0 ; 
ref.img = mrests ; 
save_untouch_nii(ref,'mrests.nii.gz') ;

refepi = load_untouch_nii('c:/shared/frefs/sumsumreg.nii.gz') ; 

for i=7:15 ; subplottight(3,3,i-6) ; 
lowrests = squeeze(mean(mrests(:,:,:,3:6),4)) ; 
plotoverlayIntensity2D(refepi.img(:,:,i),mat2gray(abs(lowrests(:,:,i))),mat2gray(lowrests(:,:,i)),270) ; 
end


%{
mmovies = squeeze(mean(avgmovies,1)) ; 
freqref.img = mmovies ; 
save_untouch_nii(freqref,'mmovies.nii.gz') ;
%}















