cd c:/shared/regute2/ ; ls 
mask = load_untouch_nii('mask.nii.gz') ; 
mask.img(:,:,1:15) = 2 ;
save_untouch_nii(mask,'mask2.nii.gz') ; 










