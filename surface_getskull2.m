clear all ; close all ; 
cd c:/shared/UTE ; ls 
s = load_untouch_nii('res_lin.nii.gz') ; 
simg = double(s.img) ; 
gs = double(simg-imfilter(simg,fspecial('gaussian',65,300))) ; 
gs2 = double(simg-imfilter(simg,fspecial('gaussian',10,5))) ; 

% k means on every slice individually
for i=1:150 ; [~,ksz(:,:,i)] = kmeans(round(mat2gray(squeeze(simg(:,:,i)))*255),2) ; end
for i=1:240 ; [~,ksx(i,:,:)] = kmeans(round(mat2gray(squeeze(simg(i,:,:)))*255),2) ; end 
for i=1:240 ; [~,ksy(:,i,:)] = kmeans(round(mat2gray(squeeze(simg(:,i,:)))*255),2) ; end
ks = ((ksy==2) .* (ksx==2) .* (ksz==2)) ;
remks = zeros(size(ks)) ;
for i=1:size(ks,3) ; 
    bw = bwconncomp(squeeze(ks(:,:,i))) ; plist = bw.PixelIdxList ; 
    maxsize = 1 ; maxind = 1 ; 
    for p=1:length(plist) ; if length(plist{p}) > maxsize ; maxind = p ; maxsize = length(plist{p}) ; end ; end 
    slice = zeros(size(ks,1),size(ks,2)) ; 
    if ~isempty(plist) ; slice(bw.PixelIdxList{maxind}) = 1 ; end
    remks(:,:,i) = slice ; 
end 
for i=1:size(ks,1) ; 
    bw = bwconncomp(squeeze(ks(i,:,:))) ; plist = bw.PixelIdxList ; 
    maxsize = 1 ; maxind = 1 ; 
    for p=1:length(plist) ; if length(plist{p}) > maxsize ; maxind = p ; maxsize = length(plist{p}) ; end ; end 
    slice = zeros(size(ks,2),size(ks,3)) ; 
    if ~isempty(plist) ; slice(bw.PixelIdxList{maxind}) = 1 ; end
    remks(i,:,:) = slice ; 
end 
 remks = medfilt3(remks) ;

dilmin = imdilate(remks,strel(ones(13,13,13)))-imdilate(remks,strel(ones(5,5,5))) ; 
bw = bwconncomp(dilmin) ; outer = zeros(size(dilmin)) ; outer(bw.PixelIdxList{1}) = 1 ; 
dil_outer = imdilate(outer,strel(ones(5,5,5))) ; 
icount = 1 ; dilprev = zeros(size(simg)) ; clear dils ; 
for i=3:2:27
   dili = imdilate(remks,strel(ones(i,i,i))) ; 
   multdil = dili.*dil_outer ; 
   if icount == 1 
      dils(icount,:,:,:) = multdil ;  
   elseif icount > 1 ;
      dils(icount,:,:,:) = (multdil - dilprev) > 0 ;      
   end
   dilprev = multdil + dilprev ; 
   icount = icount + 1 ;
end


clear vqs ; 
max_xp = 2.5 ; min_xp = -2.5 ; max_yp = 2.5 ; min_yp = -2.5 ; 
for i=1:size(dils,1) ; indimg = squeeze(dils(i,:,:,:)) ; 
inds = find(indimg==1) ;  
[x,y,z] = ind2sub(size(indimg),inds) ; 
% find the center of mass of the skull image:
[cx,cy,cz] = centmass3(simg) ; 
xdiffs = x-cx ; ydiffs = y-cy ; zdiffs = z-cz ; 
[th,phi,r] = cart2sph(xdiffs,ydiffs,zdiffs) ; 
thbrain = zeros(size(simg)) ; thbrain(inds) = th ; 
phibrain = zeros(size(simg)) ; phibrain(inds) = phi ; 
rbrain = zeros(size(simg)) ; rbrain(inds) = r ; 
[xp,yp] = pol2cart(th,max(phi)-phi) ; % theta and rho are the phi and theta brains respectively
nsteps = 550 ; xsteps = min_xp:(max_xp-min_xp)/nsteps:max_xp ; ysteps = min_yp:(max_yp-min_yp)/nsteps:max_yp ;
[xg,yg] = meshgrid(xsteps,ysteps) ; 
vq = griddata(xp,yp,double(gs2(inds)),xg,yg) ; vq(isnan(vq)) = 0 ; 
vqs(i,:,:) = vq ; 
end


for i=1:size(vqs,1) ; figure ; imagesc(squeeze(vqs(i,:,:))) ; colormap gray ; end
for i=1:3:60
    figure ; imagesc(squeeze(vqs(3,:,:))-imfilter(squeeze(vqs(3,:,:)),fspecial('gaussian',i,i))) ; 
    
end

for i=1:7 ;
    figure ;
    imagesc(squeeze(vqs(i,:,:))-imfilter(squeeze(vqs(i,:,:)),fspecial('gaussian',60,60))) ; colormap gray ;
    alldiffs(i,:,:) = squeeze(vqs(i,:,:))-imfilter(squeeze(vqs(i,:,:)),fspecial('gaussian',60,60)) ; 
end
dimg = (squeeze(mean(alldiffs(2:5,:,:)))) ; 




rgbs(:,:,1) = uint8(mat2gray(squeeze(mean(vqs(1:2,:,:)))*255) ; 
rgbs(:,:,2) = uint8(mat2gray(squeeze(mean(vqs(3:4,:,:)))*255) ;  
rgbs(:,:,3) = uint8(mat2gray(squeeze(mean(vqs(5:6,:,:)))*255) ;  








