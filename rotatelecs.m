clear all ; close all ; 
cd c:/shared/glaucoma/laurent/DICOM ; 
%cd c:/shared/lastute/alex
mask = load_untouch_nii('ute_mask.nii.gz') ; 
resute = load_untouch_nii('res_ute_1.nii.gz') ; 
l = load_untouch_nii('layers.nii.gz') ; 

maskim = pad3d(mask.img > 0,10) ;
layerim = pad3d(l.img,10) ; 
layerim(layerim==0) = nan ; 
layerim = max(max(max(layerim)))-layerim ; 
layerim(isnan(layerim)) = 0 ; 
maskshell = imdilate(maskim,strel(ones(15,15,15))) - imerode(maskim,strel(ones(9,9,9))) ; 
uteim = pad3d(resute.img,10) ; 
uteim = uteim - imfilter(uteim,fspecial('gaussian',45,35)) ; 
smoothmask = imfilter(double(maskim),fspecial('gaussian',15,5)) ; 
[gy,gx,gz] = gradient(smoothmask) ; 
rgb(:,:,:,1) = gx ;
rgb(:,:,:,2) = gy ; 
rgb(:,:,:,3) = gz ; 
rgb = (rgb) ; 
%rgb = rgb.*double(repmat(maskshell,[1,1,1,3])) ; 
locs = load('mricoords_2.mat') ; locs = locs.mricoords + 10 ; 
%{
nimg = zeros(size(mask.img)) ; 
for i=1:size(locs,2) 
    nimg(locs(1,i),locs(2,i),locs(3,i)) = 1 ; 
end
nimg = imdilate(nimg,strel(ones(5,5,5))) ; 
mask.img = nimg ; 
save_untouch_nii(mask,'locimg.nii.gz') ; 
%}

shell2 = imdilate(maskim,strel(ones(11,11,11))) - imdilate(maskim,strel(ones(5,5,5))) ; 
[gx2,gy2,gz2] = gradient(uteim) ; 
shellrgb(:,:,:,1) = gx2 ; shellrgb(:,:,:,2) = gy2 ; shellrgb(:,:,:,3) = gz2 ; 
rgbmask = repmat(shell2,[1,1,1,3]).*shellrgb ; 
elecint = shell2.*pad3d(resute.img,10) ; 

[gx,gy,gz] = meshgrid(-8:8,-8:8,-8:8) ; 
sph = double(sqrt(gx.^2 + gy.^2 + gz.^2) < 3) ; 
nz = length(find(sph==0)) ; 
sumsph = sum(sum(sum(sph))) ; 
negvoxval = -sumsph./nz ; 
sph(sph==0) = negvoxval ; 
sphute = imfilter(double(uteim),double(sph)) ; 
sphuteorig = sphute ; 
sphute = mat2gray(sphuteorig).*double(maskshell) ; %.*layerim ; 
shellute = uteim.*double(maskshell) ; 
esize = 8 ; 
for e=1:length(locs)
loc1 = locs(:,e) ; 
allrgbs(e,:,:,:,:) = squeeze(rgb(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize,:)) ; 
rgbmasked(e,:,:,:,:) = squeeze(rgbmask(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize,:)) ; 
finalint(e,:,:,:) = squeeze(elecint(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize)) ; 
tr_rgb(e,:,:,:) =  (sum(allrgbs(e,:,:,:,:),5)) ; 
rawute(e,:,:,:) = squeeze(shellute(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize)) ;
%allsphs(e,:,:,:) = squeeze(sphute(loc1(1)-10:loc1(1)+10,loc1(2)-10:loc1(2)+10,loc1(3)-10:loc1(3)+10)) ; 
end

for i=1:size(tr_rgb,1) ; 
    tri = squeeze(rawute(i,:,:,:)) ; 
    inds = find(tri~=0) ; 
    [kc,km] = kmeanscustom(uint8(mat2gray(tri(inds))*255),10) ; 
    ztri = zeros(size(tri)) ; 
    ztri(inds) = (km) ; 
    allk(i,:,:,:) = ztri ;   
    
    skullinds = find(km>=8) ; 
    [skx,sky,skz] = ind2sub(size(tri),skullinds) ; 
    clear rgbvecs
    for sk=1:length(skx)
        rgbvecs(sk,:) = squeeze(allrgbs(i,skx(sk),sky(sk),skz(sk),:)) ;        
    end   
    rgbvecs(isnan(rgbvecs)) = 0 ; 
    allvecs(i,:) = mean(rgbvecs) ; 
end

allvecs = allvecs./repmat(sqrt(sum(allvecs.^2,2)),[1,3]) ; 
for i=1:65 ; 
    figure ; for j=1:size(allk,4) ; subplot(5,5,j) ; imagesc(squeeze((allk(i,:,:,j))),[1,10]) ; end ; 
    suptitle(['x=',num2str(allvecs(i,1)),' y=',num2str(allvecs(i,2)),' z=',num2str(allvecs(i,3))]) ; 
end

