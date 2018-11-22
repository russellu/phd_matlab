%clear all ; close all 
cd c:/shared/t1atlas ; ls 
brains = dir('*in_mean_warped*gz') ; 
for b=1:length(brains)
    
    bb = load_untouch_nii(brains(b).name) ; 
    allbrains(b,:,:,:) = bb.img ; 
    
    
end

mb = squeeze(mean(allbrains,1)) ; 
bb.img = mb ; 
save_untouch_nii(bb,'sumsumreg.nii.gz') ; 


