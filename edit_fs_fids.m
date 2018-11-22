cd C:\shared\freesurfer_segs\sub1_alex\mri

t1 = load_untouch_nii('T1.mgz.nii.gz') ; 

t1.img(150,10,60) = 50 ; t1.img(154,15,110) = 150 ; t1.img(158,24,170) = 250 ; 
save_untouch_nii(t1,'T1.mgz.nii.gz') ; 

