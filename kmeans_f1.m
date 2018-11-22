cd c:/shared/claudie/sujet9 ; ls 

nii = load_untouch_nii('1mm_f1_post.nii.gz') ; 
sq = nii.img ; 
[kc,km] = kmeanscustom(uint8(mat2gray(reshape(sq,[1,numel(sq)]))*255),5) ; 
kb = reshape(km,size(nii.img)) ; 
kb(kb<5) = 0 ; 
nii.img = kb ; save_untouch_nii(nii,'kb.nii.gz') ; 
