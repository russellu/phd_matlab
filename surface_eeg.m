cd c:/shared/UTE ; ls 
s = load_untouch_nii('MONG_01_RB_DelRec_-_UTE_EEG_3_3.nii') ; 
simg = s.img ; 
%res = reshape(simg,[size(simg,2),size(simg,2)*size(simg,3)]) ; 
%eq = localhist(res,20,.5) ; reseq = reshape(eq,size(simg)) ; 
gs = simg-imfilter(simg,fspecial('gaussian',45,300)) ; 
thresh = graythresh(gs) ; 
skull = mat2gray(gs)>(thresh/1.2) ; 
skullext = (imdilate(skull,strel(ones(5,5,5)))) ; skullext2 = imdilate(skullext,strel(ones(3,3,3))) ; 
border = (skullext-skullext2) == -1 ; 
comps = bwconncomp(border)  ;
plist = comps.PixelIdxList{1} ; 
skull2 = zeros(size(skull)) ; skull2(plist) = 1 ; 

dilskull = imdilate(skull,strel(ones(3,3,3))) ; 
inner = imdilate(skull2,strel(ones(7,7,7))).*dilskull ; 

sk = imdilate(skull,strel(ones(5,5,5))).*imdilate(skull2,strel(ones(5,5,5))) ; 

inds = find(sk==1) ;  
[x,y,z] = ind2sub(size(skull),inds) ; 

% find the center of mass of the skull image:
[cx,cy,cz] = centmass3(gs) ; 
xdiffs = x-cx ; ydiffs = y-cy ; zdiffs = z-cz ; 
[th,phi,r] = cart2sph(xdiffs,ydiffs,zdiffs) ; 

thbrain = zeros(size(simg)) ; thbrain(inds) = th ; 
phibrain = zeros(size(simg)) ; phibrain(inds) = phi ; 
rbrain = zeros(size(simg)) ; rbrain(inds) = r ; 

% thbrain is the radial component, phibrain is the angular (for the pancake
% view)
[xp,yp] = pol2cart(th,max(phi)-phi) ; % theta and rho are the phi and theta brains respectively

nsteps = 400 ; xsteps = min(xp):(max(xp)-min(xp))/nsteps:max(xp) ; ysteps = min(yp):(max(yp)-min(yp))/nsteps:max(yp) ;
[xg,yg] = meshgrid(xsteps,ysteps) ; 
vq = griddata(xp,yp,double(gs(inds)),xg,yg) ; vq(isnan(vq)) = 0 ; 



% try some sharpening filters of different radii
for rad1=0:3:20 ;
    for rad2=21:3:50
%rad1 = 10 ; rad2= 30 ;
[gridx,gridy] = meshgrid(-size(vq,1)/2:size(vq,1)/2-1,-size(vq,2)/2:size(vq,2)/2-1) ; 
circle = sqrt(gridx.^2+gridy.^2)<rad2 & sqrt(gridx.^2+gridy.^2)>rad1 ; 
filtimg = abs(ifft2(fftshift(fft2(vq)).*circle)) ; 
figure,
imagesc(filtimg) ; 

    end
end





