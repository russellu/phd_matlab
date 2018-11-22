clear all ; close all ; 
sujets = {'sujet1','sujet2','sujet3','sujet4','sujet5','sujet6','sujet8','sujet9','sujet10','sujet11','sujet12'} ; 
for s=1:length(sujets) ; 
cd(['c:/shared/claudie/',sujets{s}]) ; ls 
mask = load_untouch_nii('automask.nii.gz') ; 
smoothmask = imfilter(double(mask.img),fspecial('gaussian',27,9)) ; 
[gy,gx,gz] = gradient(smoothmask) ; 
clear rgb
rgb(:,:,:,1) = gx ; rgb(:,:,:,2) = gy ; rgb(:,:,:,3) = gz ; 

parietal = load_untouch_nii('parietal_zone.nii.gz') ; 
[cx,cy,cz] = centmass3(parietal.img) ; 

clear vecs ; 
inds = find(imdilate(parietal.img,strel(ones(3,3,3)))==1) ; 
[subx,suby,subz] = ind2sub(size(parietal.img),inds) ; 
for i=1:length(inds)
    vecs(i,:) = squeeze(rgb(subx(i),suby(i),subz(i),:)) ; 
end

normvecs = vecs./repmat(sqrt(sum(vecs.^2,2)),[1,3]) ; 

mvec = mean(normvecs) ; 

voxels = repmat(1:60,[3,1])' ; 
linevals = 1:60 ; 
voxelsubs = round(voxels.*repmat(mvec,[60,1])) ; 
voxelsubs(:,1) = voxelsubs(:,1) + cx ; voxelsubs(:,2) = voxelsubs(:,2) + cy ; voxelsubs(:,3) = voxelsubs(:,3) + cz ; 

gray = load_untouch_nii('ss.nii.gz') ; gray.img = double(gray.img>0) ; 

lineinds = sub2ind(size(parietal.img),voxelsubs(:,1),voxelsubs(:,2),voxelsubs(:,3)) ;
linebrain = zeros(size(parietal.img)) ; 
ordbrain = zeros(size(parietal.img)) ; 
linebrain(lineinds) = 1 ; ordbrain(lineinds) = linevals ; 
parietal.img = imdilate(linebrain,strel(ones(3,3,3))) ; save_untouch_nii(parietal,'linebrain.nii.gz') ; 
parietal.img = ordbrain ; save_untouch_nii(parietal,'ordbrain.nii.gz') ; 

grayord = gray.img.*ordbrain ; 
grayinds = find(grayord~=0) ; 
grayindvals = ordbrain(grayinds) ; 
min_gray = min(grayindvals) ; 

sphcount = 1 ; 
uniques = unique(ordbrain) ; uniques(1) = [] ; 
for mg=1:20
    spherebrain = zeros(size(ordbrain)) ; 
    [cx,cy,cz] = ind2sub(size(ordbrain),find(ordbrain==uniques(mg))) ; 
    [gx,gy,gz] = ndgrid(-cx:size(ordbrain,1)-cx-1,-cy:size(ordbrain,2)-cy-1,-cz:size(ordbrain,3)-cz-1) ; 
    sph = sqrt(gx.^2 + gy.^2 + gz.^2) < 10 ; 
    spherebrain((sph==1)) = mg ; 
    parietal.img = spherebrain ; 
    save_untouch_nii(parietal,['spheres/sphere_',num2str(sphcount),'.nii.gz']) ; 
    sphcount = sphcount + 1  ; 
end

end

