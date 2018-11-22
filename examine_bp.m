cd c:/shared/badger_mri/alex/nii ; ls ; 
bp = load_untouch_nii('bp_reg_topup_mc_retino_gamma_01.nii.gz') ;
reg = load_untouch_nii('reg_topup_mc_retino_gamma_01.nii.gz') ; 

meanreg = squeeze(mean(reg.img,4)) ; 

vox = [32,52,6] ;