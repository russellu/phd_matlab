clear all ; close all ; 
cd c:/shared/UTE ; ls 
s = load_untouch_nii('MONG_01_RB_DelRec_-_UTE_EEG_3_3.nii') ; 
simg = s.img ; 
gs = double(simg-imfilter(simg,fspecial('gaussian',65,300))) ; 
[~,km] = kmeans(uint8(mat2gray(reshape(gs,[size(gs,1),size(gs,2)*size(gs,3)]))*255),2) ; 
kim = reshape(km,size(gs)) ; 
medk = medfilt3(kim==2) ; bw = bwconncomp(medk) ; medk = zeros(size(medk)) ; medk(bw.PixelIdxList{1}) = 1 ; 
outer_strip = imdilate(medk,strel(ones(7,7,7)))-imdilate(medk,strel(ones(3,3,3))) ; bw = bwconncomp(outer_strip) ; 
outer_strip = zeros(size(outer_strip)) ; outer_strip(bw.PixelIdxList{1}) = 1 ; 
dil_outer = imdilate(outer_strip,strel(ones(3,3,3))) ; 

icount = 1 ; dilprev = zeros(size(simg)) ; clear dils ; 
for i=3:2:9
   dili = imdilate(medk,strel(ones(i,i,i))) ; 
   multdil = dili.*dil_outer ; 
   if icount == 1 
      dils(icount,:,:,:) = multdil ;  
   elseif icount > 1 ;
      dils(icount,:,:,:) = (multdil - dilprev) > 0 ;      
   end
   dilprev = multdil + dilprev ; 
   icount = icount + 1 ;
end

sumdils = squeeze(sum(dils,1)) ; vals = gs(sumdils==1) ; 
[~,km] = kmeans(uint8(mat2gray(vals)*255),5) ; 
kdisp = zeros(size(gs)) ; kdisp(sumdils==1) = km ; 

clear vqs ; 
max_xp = 2.5 ; min_xp = -2.5 ; max_yp = 2.5 ; min_yp = -2.5 ; 
for i=1:4 ; indimg = squeeze(dils(i,:,:,:)) ; 
%indimg = squeeze(dils(1,:,:,:)+dils(2,:,:,:))>0 ;
%indimg = sim ; 
inds = find(indimg==1) ;  
[x,y,z] = ind2sub(size(indimg),inds) ; 
% find the center of mass of the skull image:
[cx,cy,cz] = centmass3(gs) ; 
xdiffs = x-cx ; ydiffs = y-cy ; zdiffs = z-cz ; 
[th,phi,r] = cart2sph(xdiffs,ydiffs,zdiffs) ; 
thbrain = zeros(size(simg)) ; thbrain(inds) = th ; 
phibrain = zeros(size(simg)) ; phibrain(inds) = phi ; 
rbrain = zeros(size(simg)) ; rbrain(inds) = r ; 
[xp,yp] = pol2cart(th,max(phi)-phi) ; % theta and rho are the phi and theta brains respectively
nsteps = 550 ; xsteps = min_xp:(max_xp-min_xp)/nsteps:max_xp ; ysteps = min_yp:(max_yp-min_yp)/nsteps:max_yp ;
[xg,yg] = meshgrid(xsteps,ysteps) ; 
vq = griddata(xp,yp,double(gs(inds)),xg,yg) ; vq(isnan(vq)) = 0 ; 
%figure,imagesc(vq) ; 
vqs(i,:,:) = vq ; 
end
for i=1:4 ; subplot(2,2,i) ; imagesc(squeeze(vqs(i,:,:))) ; colormap gray ; end

elecs = (squeeze(mean(vqs(1:2,:,:)))); 

%for i=1:4 ; s.img = squeeze(dils(i,:,:,:)) ; save_untouch_nii(s,['dil_',num2str(i),'.nii.gz']) ; end
%imagesc(squeeze(vqs(1,:,:)+(vqs(2,:,:)-vqs(4,:,:))))










