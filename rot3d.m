[x,y,z] = ndgrid(-1:.025:1) ;
%blob = z <= 0 & z >= -0.75 & x.^2 + y.^2 <= sqrt(0.25);
%blob = blob | (z > 0 & (abs(x) + abs(y) <= (0.5 - z)));
for i=1:length(allvecs) ; if isnan(allvecs(i,1)) ; allvecs(i,:) = [1,0,0] ; end ; end
for elecn = 1:65 

blob = squeeze(allk(elecn,:,:,:)) ; 
shellblob = squeeze(finalint(elecn,:,:,:)) ; 
rawuteblob = squeeze(rawute(elecn,:,:,:)) ; 

%{
figure,
p = patch(isosurface(blob,0.5));
set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');
daspect([1 1 1]);
view(3)
camlight
lighting gouraud
%}
blob_center = (size(blob) + 1) / 2 ; 

T1 = [1 0 0 0
      0 1 0 0
      0 0 1 0
    -blob_center 1] ; 

for theta=2%=-pi:.02:pi
%T2 = [cos(theta) ,0 ,-sin(theta),0 ;0 ,1,0,0;sin(theta),0,cos(theta),0; 0,0,0,1] ; 
a=  theta ; 
rotx = [1,0,0,0 ; 0,cos(a),-sin(a),0 ; 0,sin(a),cos(a),0 ; 0,0,0,1] ;
roty = [cos(a),0,sin(a),0 ; 0,1,0,0 ; -sin(a),0,cos(a),0 ; 0,0,0,1] ; 
rotz = [cos(a),-sin(a),0,0 ; sin(a),cos(a),0,0 ; 0,0,1,0 ; 0,0,0,1] ; 
T2 = rotz ; 
T5 = rotx ;  

rot = vrrotvec2mat(vrrotvec(allvecs(22,:),allvecs(elecn,:))) ;

rot(4,4) = 1 ; 

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
rawuteblob2 = tformarray(rawuteblob, tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F);

blobred = tformarray(squeeze(rgbmasked(elecn,:,:,:,1)), tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F);
blobgreen = tformarray(squeeze(rgbmasked(elecn,:,:,:,2)), tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F);
blobblue = tformarray(squeeze(rgbmasked(elecn,:,:,:,3)), tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F);

[gy,gx,gz] = gradient(rawuteblob2) ;
elecrgb(:,:,:,1) = gx ; elecrgb(:,:,:,2) = gy ; elecrgb(:,:,:,3) = gz ; 
%elecrgb = elecrgb./repmat(sum(elecrgb.^2,4),[1,1,1,3]) ; 
elecn
%figure,
%for i=1:15 ; subplot(5,5,i); imagesc(mat2gray(squeeze(elecrgb(:,:,i,:)))) ; end

%for i=1:size(blob2,3) ; subplot(5,5,i) ; imagesc((squeeze(rawuteblob2(:,:,i))),[0,1000]) ; end 


allblobs(elecn,:,:,:) = blob2 ; 
%clf
%{
p = patch(isosurface(squeeze(mean(allblobs(:,:,:,:),1)),0.5));
set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
daspect([1 1 1]);
view(3)
camlight
lighting gouraud
%}
%getframe  
end

end
%{
mblob = squeeze(mean(allblobs)) ; 

for i=1:65
   convblobs(i,:,:,:) = imfilter(squeeze(allblobs(i,:,:,:)),mblob(5:end-5,5:end-5,5:end-5)) ;  
end
for i=1:65 ; figure ; for j=1:25 ; subplot(5,5,j) ; imagesc(squeeze(newblobs(i,:,:,j)),[0,800]) ; end ; end
for i=1:65 ; [cx(i),cy(i),cz(i)] = ind2sub(size(squeeze(convblobs(1,:,:,:))),find(squeeze(convblobs(i,:,:,:))==max(max(max(convblobs(i,:,:,:)))))) ; end
for i=1:size(allblobs,1) ; newblobs(i,:,:,:) = allblobs(i,:,:,:) ; newblobs(i,cx(i)-1:cx(i)+1,cy(i)-1:cy(i)+1,cz(i)-1:cz(i)+1) = 600 ; end ;  

[circx,circy] = meshgrid(-size(allblobs,2)/2:size(allblobs,2)/2-1,-size(allblobs,3)/2:size(allblobs,3)/2-1) ; 
circmask = sqrt(circx.^2 + circy.^2) < 3 ; 
sumcirc = sum(sum(circmask)) ; 
circmask = circmask./sumcirc ; 
sumnotcirc = sum(sum(circmask==0)) ; 
negval = -1/sumnotcirc ; 
circmask(circmask==0) = negval ; 

for i=1:size(allblobs,1)
   for j=1:size(allblobs,4)
       convblobs(i,:,:,j) = imfilter(squeeze(allblobs(i,:,:,j)),circmask) ;      
   end
end
%}

for i=1:65
blobi = squeeze(allblobs(i,:,:,:)) ; 
red = squeeze(rgbmasked(i,:,:,:,1)) ; green = squeeze(rgbmasked(i,:,:,:,2)) ; blue = squeeze(rgbmasked(i,:,:,:,3)) ; 
blobinds = find(blobi~=0) ; 
[bx,by,bz] = ndgrid(-size(blobi,1)/2:size(blobi,1)/2-1,-size(blobi,2)/2:size(blobi,2)/2-1,-size(blobi,3)/2:size(blobi,3)/2-1) ; 
[th,phi,r] = cart2sph(bx,by,bz) ; 
blobvals = blobi(blobinds) ; 
blobx = th(blobinds) ; bloby = phi(blobinds) ; blobz = r(blobinds) ; 
blobx2 = bx(blobinds) ; bloby2 = by(blobinds) ; blobz2 = bz(blobinds) ; 

blobvec = [mat2gray(blobvals)*4,mat2gray(blobx),mat2gray(bloby),mat2gray(blobz)] ; 
km = kmeans(blobvec,4) ; 
newblob = zeros(size(blobi)) ; 
newblob(blobinds) = km ; 
newblobs(i,:,:,:) = medfilt3(newblob) ; 

% examine the clusters for electrode-ness.
end
for i=1:65 ; figure ; for j=1:25 ; subplot(5,5,j) ; imagesc(squeeze(newblobs(i,:,:,j)),[0,]) ; end ; end
for i=1:65 ; subplot(6,11,i) ; imagesc(squeeze(max(medfilt3(squeeze(allblobs(i,:,:,:))),[],1))) ; end



[gx,gy,gz] = ndgrid(-25:25,-25:25,-25:25) ; 
[th,phi,r] = cart2sph(gx,gy,gz) ; 
shell = r>9 & r<10 ; 
inds = find(shell==1) ; 
vecs  = [gx(inds),gy(inds),gz(inds)] ; 
vecs = vecs./repmat(sum(vecs.^2,2),[1,3]) ; 

for ec=1:65 ; figure,
    vecblobs = zeros(size(vecs,1),size(finalint,2),size(finalint,3),size(finalint,4)) ; 
    for v=1:length(vecs) ; 
        shellblob = squeeze(finalint(ec,:,:,:)) ; 
        rot = vrrotvec2mat(vrrotvec(vecs(v,:),allvecs(ec,:))) ;
        rot(4,4) = 1 ; 
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
        vecblobs(v,:,:,:) = blob2 ; v,
    end
    maxvecblobs = max(vecblobs,[],4) ; 
    icount =1 ; for i=1:10:380 ; subplot(6,7,icount) ; imagesc(squeeze(maxvecblobs(i,:,:))) ; icount = icount + 1 ; end 

    for i=1:size(maxvecblobs,1) ; repblobs(i,:,:,:) = repmat(squeeze(maxvecblobs(i,:,:)),[1,1,25]) ; end
    invblobs = zeros(size(vecs,1),size(finalint,2),size(finalint,3),size(finalint,4)) ; 
    for i=1:size(vecs,1)
        shellblob = squeeze(repblobs(i,:,:,:)) ; 
        rot = vrrotvec2mat(vrrotvec(vecs(i,:),allvecs(ec,:))) ;
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

for i=1:65 ; figure
    for j=1:17 ; subplot(5,5,j) ; imagesc(squeeze(allinvs(i,:,:,j)),[0,800]) ; end
end

clear eleck
for i=1:65 ; 
    ki = (squeeze(allinvs(i,:,:,:))) ; 
    %[kc,km] = kmeanscustom(uint8(mat2gray(reshape(ki,[1,numel(ki)]))*255),10) ; 
    %eleck(i,:,:,:) = reshape(km,size(ki)) ;  
    [tv,topinds] = sort(ki(:),'descend') ; 
    newki = zeros(size(ki)) ; 
    newki(topinds(1:100)) = 1 ; 
    newki2 = zeros(size(ki)) ; 
    newki2(topinds(101:200)) = 1 ; 
    eleck(i,:,:,:) = newki ;
    eleck2(i,:,:,:) = newki2 ; 
end


icount = 1 ; 
for i=1:size(eleck,1)
    ki = squeeze(eleck(i,:,:,:)>9) ; 
    bw = bwconncomp(ki) ; 
    pix = bw.PixelIdxList ; 
    if length(pix) > 1
        icount = icount +1  ;        
    end
end
icount/65 



elecimg = zeros(size(shellute)) ; 
segimg2 = zeros(size(shellute)) ; 
segimg = zeros(size(shellute)) ; 
boximg = zeros(size(shellute)) ; 
coordimg = zeros(size(shellute)) ; 
esize = 8 ; 
for e=1:length(locs)
loc1 = locs(:,e) ; 
ek = squeeze(eleck(e,:,:,:)) ; ek2 = squeeze(eleck2(e,:,:,:)) ; 

[cx,cy,cz] = centmass3(ek) ; 
coord = zeros(size(ek)) ; coord(cx,cy,cz) = e ; 
zimg = zeros(size(ek)) ; zimg(cx,cy,cz) = 1 ; zimg = imdilate(zimg,strel(ones(3,3,3))) ; 
elecimg(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize) = zimg ; 
segimg(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize) = ek ; 
segimg2(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize) = ek2 ; 

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











%{
elecimg = zeros(size(shellute)) ; 
for i=1:65 ; 
    loc1 = locs(:,i) ; 
   inti = squeeze(finalint(i,:,:,:)) ; 
   inds = find(inti~=0) ; 
   vals = inti(inds) ; 
   [kc,km] = kmeanscustom(uint8(mat2gray(vals)*255),5) ; 
   newint = zeros(size(inti)) ; 
   newint(inds) = km ; 
   allks(i,:,:,:) = newint ; 
   elecimg(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize) = newint>2 & newint<4 ; 

end
%for i=1:65 ; figure ; for j=1:17 ; subplot(4,5,j) ; imagesc(squeeze(allks(i,:,:,j)),[0,5]) ; end ; end
resute.img = uint8(elecimg(10:end-11,10:end-11,10:end-11)) ; 
save_untouch_nii(resute,'resute_pts.nii.gz') ; 
%}










