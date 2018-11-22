clear all ; close all ; 
cd c:/shared/lastute/alex ; 
locs = load('mricoords_1.mat') ; locs = locs.mricoords ; 

ute = load_untouch_nii('res_ute.nii.gz') ; 
mask = load_untouch_nii('fnirt/ute_mask.nii.gz') ; 
allnorms = zeros(size(ute)) ; 

for s=1:length(locs) ; 
loc1 = locs(:,s) ;
pimg = zeros(size(ute.img)) ; 

pimg(loc1(1),loc1(2),loc1(3)) = 1 ; 
inds = find(imdilate(pimg,strel(ones(12,12,12))).*mask.img) ; 

[cx,cy,cz] = ind2sub(size(ute.img),inds) ; 
[gx,gy,gz] = gradient(double(ute.img)) ; 

rgb(:,:,:,1) = (-gx).*mask.img ; 
rgb(:,:,:,2) = (-gy).*mask.img ; 
rgb(:,:,:,3) = (-gz).*mask.img ; 

for i=1:length(inds)
   vecs(i,:) = rgb(cx(i),cy(i),cz(i),:) ;  
end

meanvec = squeeze(mean(vecs)) ; 
normvec = meanvec./sqrt(sum(meanvec.^2)) ; 
[mx,my,mz] = centmass3(pimg) ; 
%[gridx,gridy,gridz] = ndgrid(-mx:size(p.img,1)-mx-1,-my:size(p.img,2)-my-1,-mz:size(p.img,3)-mz-1) ; 
v = repmat(1:50,[3,1]) ; 
repnorm = repmat(normvec,[50,1])' ; 
vecvals = repnorm.*v ; 
vecvals(1,:) = vecvals(1,:) + mx ; vecvals(2,:) = vecvals(2,:) + my ; vecvals(3,:) = vecvals(3,:) + mz ; 
vecvals = round(vecvals) ; 
normbrain = zeros(size(pimg)) ; 
for i=1:size(vecvals,2)
    
    normbrain(vecvals(1,i),vecvals(2,i),vecvals(3,i)) = 1 ;  
end


allnorms = allnorms + normbrain(1:240,1:240,1:170) ; 
end
mask.img = allnorms ; save_untouch_nii(mask,'normbrain.nii.gz') ; 




