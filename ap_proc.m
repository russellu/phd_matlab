clear all ; close all
cd c:/shared/highres3d/ ; ls 
cat = load_untouch_nii('cat.nii') ; 

img = (mean(mean(cat.img,3),4)>0) ; 
goods = find(sum(img,2)>0) ; 

newimg = cat.img(goods,:,:,:) ; 
save_nii(make_nii(newimg),'newimg.nii.gz') ; 


cd c:/shared/highres3d/ ; ls 

swi = load_untouch_nii('20170504_124948WIPSWIpSENSEs1601a1016.nii.gz') ; 
magswi = zeros(size(swi.img,1),size(swi.img,2),size(swi.img,3)) ; 
magswi(:,:,1:2:end) = swi.img(:,:,1:108,1) ; 
magswi(:,:,2:2:end) = swi.img(:,:,1:108,2) ;
phaseswi = zeros(size(swi.img,1),size(swi.img,2),size(swi.img,3)) ; 
phaseswi(:,:,1:2:end) = swi.img(:,:,109:end,1) ; 
phaseswi(:,:,2:2:end) = swi.img(:,:,109:end,2) ;
save_nii(make_nii(magswi),'magswi.nii.gz') ; 
save_nii(make_nii(phaseswi),'phaseswi.nii.gz') ; 

save_nii(make_nii(phaseswi.*magswi),'bothswi.nii.gz') ; 


cd mel ; 
meanf = load_untouch_nii('mean.nii.gz') ; 
meanf.img = meanf.img - imfilter(meanf.img,fspecial('gaussian',35,35)) ; 
save_untouch_nii(meanf,'meanf.nii.gz') ; 

cd ../
swi = load_untouch_nii('swi_in_t1.nii.gz') ; 
newswi = zeros(size(swi.img)) ; 
newswi(100:200,1:250,50:130) = swi.img(100:200,1:250,50:130) ; 
swi.img = newswi ; save_untouch_nii(swi,'cropswi.nii.gz') ; 



