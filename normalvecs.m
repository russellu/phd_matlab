cd c:/shared/ATLASES ; 

p = load_untouch_nii('parietal.nii.gz') ; 
mni = load_untouch_nii('MNI152_T1_1mm.nii.gz') ; 

inds = find(imdilate(p.img,strel(ones(7,7,7)))) ; 
[cx,cy,cz] = ind2sub(size(p.img),inds) ; 
[gx,gy,gz] = gradient(double(mni.img)) ; 
rgb(:,:,:,1) = (gx) ; 
rgb(:,:,:,2) = (gy) ; 
rgb(:,:,:,3) = (gz) ; 

for i=1:length(inds)
   vecs(i,:) = rgb(cx(i),cy(i),cz(i),:) ;  
end

meanvec = squeeze(mean(vecs)) ; 
normvec = meanvec./sqrt(sum(meanvec.^2)) ; 
[mx,my,mz] = centmass3(p.img) ; 
%[gridx,gridy,gridz] = ndgrid(-mx:size(p.img,1)-mx-1,-my:size(p.img,2)-my-1,-mz:size(p.img,3)-mz-1) ; 
v = repmat(1:100,[3,1]) ; 
repnorm = repmat(normvec,[100,1])' ; 
vecvals = repnorm.*v ; 
vecvals(1,:) = vecvals(1,:) + mx ; vecvals(2,:) = vecvals(2,:) + my ; vecvals(3,:) = vecvals(3,:) + mz ; 
vecvals = round(vecvals) ; 
normbrain = zeros(size(p.img)) ; 
for i=1:size(vecvals,2)
    normbrain(vecvals(1,i),vecvals(2,i),vecvals(3,i)) = 1 ;  
end

p.img = imdilate(normbrain,strel(ones(3,3,3))) ; save_untouch_nii(p,'normbrain.nii.gz') ; 





