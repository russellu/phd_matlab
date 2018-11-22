clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 

for sb=1:length(subs)
    
    
    cd(['E:\rawbadger\badger_mri\',subs{sb},'\nii\warped']);   
    corrs = dir('eegcorrs*'); 
    for corr=1:length(corrs)
       nii = load_untouch_nii(corrs(corr).name); 
       allniis(:,:,:,:,sb,corr) = nii.img; 
        
        
    end
            
end


mean_all = squeeze(mean(mean(allniis,5),6)); 
mean_rest = squeeze(mean(allniis(:,:,:,:,:,6),5));
mean_movie = squeeze(mean(allniis(:,:,:,:,:,5),5)); 
mean_gamma = squeeze(mean(allniis(:,:,:,:,:,3:4),5)); 
mean_retino = squeeze(mean(allniis(:,:,:,:,:,1:2),5)); 


cd E:\meanepis
nii.img = mean_all; save_untouch_nii(nii,'mean_hrf_all.nii.gz'); 
nii.img = mean_rest; save_untouch_nii(nii,'mean_hrf_rest.nii.gz'); 
nii.img = mean_movie; save_untouch_nii(nii,'mean_hrf_movie.nii.gz'); 
nii.img = mean_gamma; save_untouch_nii(nii,'mean_hrf_gamma.nii.gz'); 
nii.img = mean_retino; save_untouch_nii(nii,'mean_hrf_retino.nii.gz'); 










