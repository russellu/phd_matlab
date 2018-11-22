cd c:/shared/testorient ; 
prefwd = load_untouch_nii('steve_pre_01.nii.gz') ;
prerev = load_untouch_nii('steve_pre_rev_01.nii.gz') ; 

cd c:/shared/steve ; 
postfwd = load_untouch_nii('steve_post_01.nii.gz') ;
postrev = load_untouch_nii('steve_post_rev_01.nii.gz') ; 

postfwd.img = permute(prefwd.img,[1,3,2,4]) ; 
postfwd.img = postfwd.img(:,:,end:-1:1,:) ; 
save_untouch_nii(postfwd,'steve_pre_01.nii.gz') ; 
postrev.img = permute(prerev.img,[1,3,2,4]) ; 
postrev.img = postrev.img(:,:,end:-1:1,:) ; 
save_untouch_nii(postrev,'steve_pre_rev_01.nii.gz') ; 