
%%% pre process for the ANTS registration
%clear all ; close all ; disp('loading data...') ; 
%cd('c:/shared/badger/felix/ute') ; s = load_untouch_nii('resamp_ute.nii.gz') ; 
%simg = s.img ; gs = double((simg)-imfilter((simg),fspecial('gaussian',90,90))) ; 
%s.img = gs ; save_untouch_nii(s,'gs.nii.gz') ; 

clear all ; close all ; disp('loading data...') ; 
cd('c:/shared/badger/felix/ute/fnirty') ; s = load_untouch_nii('gs.nii.gz') ; 
elayers = load_untouch_nii('warp_layers.nii.gz') ;  eimg = elayers.img ;   
for i=1:8 ; dils(i,:,:,:) = (eimg==i) ; end

s = load_untouch_nii('gs.nii.gz') ; 
simg = s.img ; %gs = double((simg)-imfilter((simg),fspecial('gaussian',15,15))) ;
gs = double(s.img) ; 
%s.img = gs ; save_untouch_nii(s,'gs.nii.gz') ; 
clear vqs ; 
max_xp = 2.5 ; min_xp = -2.5 ; max_yp = 2.5 ; min_yp = -2.5 ; 
for i=1:size(dils,1) ; indimg = squeeze(dils(i,:,:,:)) ; 
inds = find(indimg==1) ;  
[x,y,z] = ind2sub(size(indimg),inds) ; 
% find the center of mass of the skull image:
[cx,cy,cz] = centmass3(simg) ; 
xdiffs = x-cx ; ydiffs = y-cy ; zdiffs = z-cz ; 
[th,phi,r] = cart2sph(xdiffs,ydiffs,zdiffs) ; 
allth{i} = th ; allphi{i} = phi ; allr{i} = r ; 
thbrain = zeros(size(simg)) ; thbrain(inds) = th ; 
phibrain = zeros(size(simg)) ; phibrain(inds) = phi ; 
rbrain = zeros(size(simg)) ; rbrain(inds) = r ; 
[xp,yp] = pol2cart(th,max(phi)-phi) ; % theta and rho are the phi and theta brains respectively
nsteps = 550 ; xsteps = min_xp:(max_xp-min_xp)/nsteps:max_xp ; ysteps = min_yp:(max_yp-min_yp)/nsteps:max_yp ;
[xg,yg] = meshgrid(xsteps,ysteps) ; 
vq = griddata(double(xp),double(yp),gs(inds),double(xg),double(yg)) ; vq(isnan(vq)) = 0 ; 
vqs(i,:,:) = vq ; 
end

for i=1:8 ; 
    dqs(i,:,:) = (squeeze(vqs(i,:,:)))-imfilter((squeeze(vqs(i,:,:))),fspecial('gaussian',20,20)) ;  
end
rgbs(:,:,3) = uint8(mat2gray(squeeze(mean(dqs(5:6,:,:),1)))*255) ;  
rgbs(:,:,2) = uint8(mat2gray(squeeze(mean(dqs(3:4,:,:),1)))*255) ;  
rgbs(:,:,1) = uint8(mat2gray(squeeze(mean(dqs(1:2,:,:),1)))*255) ;  
figure,imagesc(rgbs)

allrgbs(el,:,:,:) = rgbs ; 
allvqs(el,:,:,:) = vqs ; 
imwrite(rgbs,['c:/shared/UTE_pngs/','felix','.png']) ;

% do stuff in java...
