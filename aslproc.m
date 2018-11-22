cd c:/shared/DICOM ; 
catasl = load_untouch_nii('catasl.nii.gz')  ;
cat1 = catasl.img(:,:,:,1:156) ; cat2 = catasl.img(:,:,:,157:end) ; 


for i=1:size(ctrl,4)
    
   ctrl = cat2(:,:,:,1:156/2) ; tag = cat2(:,:,:,158/2:end) ; 
   inds = 1:size(ctrl,4) ; 
   taginds = find(inds-i >=-2 & inds-i <=2) ; 
   meani = mean(tag1(:,:,:,taginds),4) ; 
   subbed1(:,:,:,i) = ctrl(:,:,:,i) - meani ; 
   
      ctrl = cat2(:,:,:,1:156/2) ; tag = cat2(:,:,:,158/2:end) ; 
   inds = 1:size(ctrl,4) ; 
   taginds = find(inds-i >=-2 & inds-i <=2) ; 
   meani = mean(tag1(:,:,:,taginds),4) ; 
   subbed1(:,:,:,i) = ctrl(:,:,:,i) - meani ; 
end




ref = load_untouch_nii('refasl.nii.gz') ; 
ref.img = (subbed1+subbed2)/2 ; 
save_untouch_nii(ref,'subbed_1.nii.gz') ; 


