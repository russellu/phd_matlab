clear all ; close all ; 
cd('C:\shared\badger\alex\Russell BADGER 2015-11-12') ; ls 
ref1 = load_untouch_nii('reg_mc_Test_Russell_2015_11_12_WIP_EEG-fMRI_MB3_3.75mm_SENSE_6_1.nii.gz') ; 
topup = load_untouch_nii('Test_Russell_2015_11_12_WIP_revEEG-fMRI_MB3_3.75mm_SENSE_12_1.nii') ; 
reptopup = repmat(topup.img,[1,1,1,size(ref1.img,4)]) ; 
ref1.img = reptopup ; save_untouch_nii(ref1,'reptopup.nii.gz') ; 



