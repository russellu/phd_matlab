clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','tegan','valerie'}; 
scans = {'*allstim*01*','*allstim*02*','*gamma*01*','*gamma*02*','*movie*','*rest*'};
freqs = {'delta','theta','alpha','beta','mid_gamma','high_gamma'};
for sb=1:length(subs) 
    cd(['e:/sim/fmri/',subs{sb}]);
    meanimg = load_untouch_nii('mean.nii.gz'); 
    mask = meanimg.img > mean(meanimg.img(:)*1.2); 
    maskvox = find(mask == 1); 
    allzsums = load('allzsums'); allzsums = allzsums.allzsums; 
    for sc=1:length(scans)

        for fr=6%:length(freqs)
            eeg = dir(['atlas_',scans{sc},freqs{fr},'*nii']);
            eeg = load_untouch_nii(eeg.name); 
            reseeg = reshape(eeg.img,[numel(eeg.img(:,:,:,1)),size(eeg.img,4)]); 
            reseegvox = reseeg(maskvox,:); 
            
            %eeg_ts = eegfiltfft(smooth(mean(reseegvox,1))',1/0.693,0.01,1); 
            %motion_ts = eegfiltfft(smooth(allzsums{sc})',1/0.693,0.01,1); 
            eeg_ts = smooth(mean(reseegvox,1)); 
            
            %corrs(sb,sc,1) = corr2(eeg_ts(30:end-30),motion_ts(30:length(eeg_ts)-30)); 
            
            
            %mot = (zscore(mean(reseegvox,1))); 
            [sv,si] = sort(eeg_ts,'descend'); 
            allmots{sc} = eeg_ts; 
            
        end   
    end
    save('allmots','allmots'); 
end




