clear all; close all; 
cd C:\shared\coupling\alex; 
fmri = load_untouch_nii('bp_clean_retino_rest.nii.gz'); 
cd den_retino_rest; 
hz1 = load_untouch_nii(['hz_',num2str(1),'.nii.gz']); 
allhzs = zeros(size(hz1.img,1),size(hz1.img,2),size(hz1.img,3),size(hz1.img,4),45); 
nums = 1:2:89; 
for i=1:length(nums)
    hzi = load_untouch_nii(['hz_',num2str(nums(i)),'.nii.gz']); 
    allhzs(:,:,:,:,i) = hzi.img; 
end
cd .. 
grid = load_untouch_nii('grid.nii.gz'); 
brainpts = find(grid.img==1); 
[px,py,pz] = ind2sub(size(grid.img),brainpts); 



