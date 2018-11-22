cd c:/shared/glaucoma/laurent/DICOM ; ls ;

nii = load_untouch_nii('res_ute_2.nii.gz') ; 
nii.img = nii.img - imfilter(nii.img,fspecial('gaussian',20,10)) ; 
save_untouch_nii(nii,'gs_20_2.nii.gz') ; 





