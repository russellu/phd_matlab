cd('c:/shared/badger/alex/Russell BADGER 2015-11-12/topup') ;
fs=dir('unwarped_reg*') ; 
for f=1%:length(fs) ; 
   fi = load_untouch_nii(fs(f).name) ; 
   mfi = squeeze(mean(fi.img,4)./std(fi.img,0,4)) ; 
   cd .. ; f1 = load_untouch_nii('f_mc_Test_Russell_2015_11_12_WIP_EEG-fMRI_MB3_3.75mm_SENSE_10_1.nii.gz') ; 
   f1.img = mfi ; save_untouch_nii(f1,'mfi.nii.gz') ; 
   %[~,km] = kmeans(uint8(mat2gray(reshape(fi.img,[1,numel(fi.img)]))*255),2) ; 
   %kimg = reshape(km,size(fi.img)) ; 
    
   
end


