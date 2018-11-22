clear all ; close all ; 
cd c:/shared/lastute/valerie
mask = load_untouch_nii('fnirt/ute_mask.nii.gz') ; 
resute = load_untouch_nii('res_ute.nii.gz') ; 
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

locs = load('mricoords_1.mat') ; locs = locs.mricoords + 10 ; 

shell2 = imdilate(maskim,strel(ones(11,11,11))) - imdilate(maskim,strel(ones(5,5,5))) ; 
[gx2,gy2,gz2] = gradient(uteim) ; 
shellrgb(:,:,:,1) = gx2 ; shellrgb(:,:,:,2) = gy2 ; shellrgb(:,:,:,3) = gz2 ; 
rgbmask = repmat(shell2,[1,1,1,3]).*shellrgb ; 
elecint = shell2.*pad3d(resute.img,10) ; 

shellute = uteim.*double(maskshell) ; 
esize = 10 ; 
clear finalint rawute
for e=1:length(locs)
loc1 = locs(:,e) ; 
finalint(e,:,:,:) = squeeze(elecint(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize)) ; 
rawute(e,:,:,:) = squeeze(shellute(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize)) ;
end




% do the denosing

[gx,gy,gz] = ndgrid(-25:25,-25:25,-25:25) ; 
[th,phi,r] = cart2sph(gx,gy,gz) ; 
shell = r>9 & r<10 ; 
inds = find(shell==1) ; 
vecs  = [gx(inds),gy(inds),gz(inds)] ; 
vecs = vecs./repmat(sum(vecs.^2,2),[1,3]) ; 
blob = squeeze(finalint(1,:,:,:)) ; 
blob_center = (size(squeeze(finalint(1,:,:,:))) + 1) / 2 ; 

clear allinvs
for ec=1:65 ; % figure,
    vecblobs = zeros(size(vecs,1),size(finalint,2),size(finalint,3),size(finalint,4)) ; 
    for v=1:length(vecs) ; 
        shellblob = squeeze(finalint(ec,:,:,:)) ; 
        rot = vrrotvec2mat(vrrotvec(vecs(v,:),[1,0,0])) ;
        rot(4,4) = 1 ; 
        T3 = [1 0 0 0
              0 1 0 0
              0 0 1 0
              blob_center 1] ; 
        T1 = [1 0 0 0
              0 1 0 0
              0 0 1 0
             -blob_center 1] ; 
        T = T1 * rot * T3 ; 
        tform = maketform('affine', T) ;
        tformfwd(blob_center, tform) ; 
        R = makeresampler('linear', 'fill') ; 
        TDIMS_A = [1 2 3];
        TDIMS_B = [1 2 3];
        TSIZE_B = size(blob);
        TMAP_B = [];
        F = 0;
        blob2 = tformarray(shellblob, tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F);
        vecblobs(v,:,:,:) = blob2 ; v,
    end
    maxvecblobs = max(vecblobs,[],4) ; 
    %icount =1 ; for i=1:10:380 ; subplot(6,7,icount) ; imagesc(squeeze(maxvecblobs(i,:,:))) ; icount = icount + 1 ; end 

    for i=1:size(maxvecblobs,1) ; repblobs(i,:,:,:) = repmat(squeeze(maxvecblobs(i,:,:)),[1,1,25]) ; end
    invblobs = zeros(size(vecs,1),size(finalint,2),size(finalint,3),size(finalint,4)) ; 
    for i=1:size(vecs,1)
        shellblob = squeeze(repblobs(i,:,:,:)) ; 
        rot = vrrotvec2mat(vrrotvec(vecs(i,:),[1,0,0])) ;
        rot(4,4) = 1 ; 
        rot = inv(rot) ; 
        T3 = [1 0 0 0
            0 1 0 0
            0 0 1 0
            blob_center 1] ; 
        T = T1 * rot * T3 ; 
        tform = maketform('affine', T) ;
        tformfwd(blob_center, tform) ; 
        R = makeresampler('linear', 'fill') ; 
        TDIMS_A = [1 2 3];
        TDIMS_B = [1 2 3];
        TSIZE_B = size(blob);
        TMAP_B = [];
        F = 0;
        blob2 = tformarray(shellblob, tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F);
        invblobs(i,:,:,:) = blob2 ; i,
    end
    allinvs(ec,:,:,:) = squeeze(mean(invblobs,1)) ; 
end

%for i=1:65 ; figure
%    for j=1:17 ; subplot(5,5,j) ; imagesc(squeeze(allinvs(i,:,:,j)),[0,800]) ; end
%end

clear eleck eleck2 dilints
for i=1:65 ; 
    ki = (squeeze(allinvs(i,:,:,:))) ; 
    [tv,topinds] = sort(ki(:),'descend') ; 
    newki = zeros(size(ki)) ; 
    newki(topinds(1:200)) = 1 ; 
    newki2 = zeros(size(ki)) ; 
    newki2(topinds(101:200)) = 1 ; 
    eleck(i,:,:,:) = newki ;
    eleck2(i,:,:,:) = newki2 ;   
    dilint = imdilate(newki,strel(ones(3,3,3))).*squeeze(rawute(i,:,:,:)) ; 
    inds = find(dilint~=0) ; [kc,km] = kmeanscustom(uint8(mat2gray(dilint(inds))*255),2) ; 
    dilint(inds) = km ; dilint = dilint==2 ; 
    dilints(i,:,:,:) = dilint ; 
end
%for i=1:65 ; figure ; for j=1:17 ; subplot(4,5,j) ; imagesc(squeeze(dilints(i,:,:,j))) ; end ; end

elecimg = zeros(size(shellute)) ; 
segimg2 = zeros(size(shellute)) ; 
segimg = zeros(size(shellute)) ; 
boximg = zeros(size(shellute)) ; 
coordimg = zeros(size(shellute)) ; 
dilimg = zeros(size(shellute)) ; 
for e=1:length(locs)
loc1 = locs(:,e) ; 
ek = squeeze(eleck(e,:,:,:)) ; ek2 = squeeze(eleck2(e,:,:,:)) ; 
[cx,cy,cz] = centmass3(ek) ; 
coord = zeros(size(ek)) ; coord(cx,cy,cz) = e ; 
zimg = zeros(size(ek)) ; zimg(cx,cy,cz) = 1 ; zimg = imdilate(zimg,strel(ones(3,3,3))) ; 
elecimg(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize) = zimg ; 
segimg(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize) = ek ; 
segimg2(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize) = ek2 ; 
dilimg(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize) = squeeze(dilints(e,:,:,:)) ; 

boximg(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize) = 1 ; 
coordimg(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize) = coord ; 
end

resute.img = elecimg(10:end-11,10:end-11,10:end-11) ; 
save_untouch_nii(resute,'resute_pts.nii.gz') ; 
resute.img = segimg(10:end-11,10:end-11,10:end-11) ; 
resute.img = uint8(resute.img) ; 
save_untouch_nii(resute,'resute_seg.nii.gz') ; 
resute.img = segimg2(10:end-11,10:end-11,10:end-11) ; 
resute.img = uint8(resute.img) ; 
save_untouch_nii(resute,'resute_seg2.nii.gz') ; 
resute.img = boximg(10:end-11,10:end-11,10:end-11) ; 
save_untouch_nii(resute,'boximg.nii.gz') ; 
resute.img = coordimg(10:end-11,10:end-11,10:end-11) ; 
save_untouch_nii(resute,'coordimg.nii.gz') ; 
resute.img = dilimg(10:end-11,10:end-11,10:end-11) ; 
save_untouch_nii(resute,'dilimg.nii.gz') ; 









