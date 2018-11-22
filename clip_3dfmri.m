cd c:/shared/test_sat ;
ls

ims = [11,14,18,20,21,25,27] ; 
for i=1:length(ims)
img = load_untouch_nii(['IM_00',num2str(ims(i)),'.nii.gz']) ; 
meanimg = mean(img.img,4) ; 
goods = find(sum(sum(meanimg,1),3)>0) ; 
newimg = img.img(:,goods,:,:) ; 
save_nii(make_nii(newimg),['cropped_',num2str(ims(i)),'.nii.gz']) ; 
end


