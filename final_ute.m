%%% process and hand label UTE images.
clear all ; close all ; 
name = 'felix' ; 
cd(['C:\shared\utef\',name]) ; ls ; 
disp('loading raw data...') ; 
rute = load_untouch_nii(['dof12_',name,'.nii.gz']) ; ruteorig = rute.img ; 
rutime = rute.img ; 

disp('filling head with ray tracing...') ; 

utemask = load_untouch_nii('c:\shared\regute\maskinv.nii.gz') ; utemask = utemask.img ; meddilate = double(utemask) ; 

l1 = imdilate(imdilate(utemask,strel(ones(3,3,3)))-utemask,strel(ones(3,3,3))) ;
vals = ruteorig(find(l1==1)) ; valinds = find(l1==1) ;
zdil = zeros(size(meddilate)) ; zdil(valinds(find(zscore(vals)>1.5))) = 1 ; 
meddilate = meddilate + zdil ; %meddilate = medfilt3(meddilate) ; 
zsub = zeros(size(meddilate)) ; zsub(valinds(find(zscore(vals)<0))) = 1 ; 
meddilate = meddilate - zsub ; meddilate = double(medfilt3(meddilate>0)) ; 
disp('creating scalp layers...') ; 
outside = bwconncomp(meddilate==0) ; 
outim = zeros(size(meddilate)) ; 
inds = outside.PixelIdxList{1} ; 
outim(inds) = 1 ; outim = double(outim) ; 
intensitylayers = zeros(size(meddilate)) ; 
prevdil = meddilate ; 
disp('performing iterative dilation...') ; 
for i=1:12 ; 
    dilmaski = imdilate(prevdil,strel(ones(3,3,3))) .* outim ; 
    clayer = (dilmaski - prevdil) > 0 ; 
    intensitylayers(clayer==1) = i ; 
    if i==1
        prevdil = dilmaski ;
    else
        prevdil = (prevdil + dilmaski) > 0 ; 
    end
end
disp('computing final layers...') ; 
layers = intensitylayers.*(meddilate==0) ; 
layers(:,:,1:18) = 0 ; 

disp('performing pancake projection...') ; 
gsimg = double(ruteorig-imfilter(rute.img,fspecial('gaussian',11,7))) ;
% hard-coded grid variables
max_xp = 2.5 ; min_xp = -2.5 ; max_yp = 2.5 ; min_yp = -2.5 ; 

for layer_index = 1:12
    layerinds = find(layers==layer_index) ; 
    [lx,ly,lz] = ind2sub(size(layers),layerinds) ; 
    [cx,cy,cz] = centmass3(layers) ; 
    xlayer_diffs = lx-cx ; ylayer_diffs = ly-cy ; zlayer_diffs = lz-cz ; 
    [theta,phi,rho] = cart2sph(xlayer_diffs,ylayer_diffs,zlayer_diffs) ; 
    if layer_index==1 ; ftheta = theta ; fphi = phi ; frho = rho ; end
    theta_brain = zeros(size(layers)) ; phi_brain = zeros(size(layers)) ; rho_brain = zeros(size(layers)) ;
    theta_brain(layerinds) = theta ; phi_brain(layerinds) = phi ; rho_brain(layerinds) = rho ; 
    [pancake_x,pancake_y] = pol2cart(theta,max(phi)-phi) ; % theta and phi are the rotation point and height point (z) of the head
    nsteps = 550 ; xsteps = min_xp:(max_xp-min_xp)/nsteps:(max_xp) ; ysteps = min_yp:(max_yp-min_yp)/nsteps:(max_yp) ;
    [xg,yg] = meshgrid(xsteps,ysteps) ; 
    vq = griddata(double(pancake_x),double(pancake_y),gsimg(layerinds),double(xg),double(yg)) ; vq(isnan(vq)) = 0 ; 
    vqs(layer_index,:,:) = vq ; 
end

elecorder = {'FP1','FPZ','FP2','AF8','AF4','GND','AF3','AF7','F7','F5','F3','F1','FZ','F2','F4','F6','F8','FT10','FT8','FC6','FC4','FC2','REF','FC1','FC3','FC5','FT7','FT9',...
    'T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP10','TP8','CP6','CP4','CP2','CPZ','CP1','CP3','CP5','TP7','TP9','P7','P5','P3','P1','PZ','P2','P4','P6','P8',...
    'PO8','PO4','POZ','PO3','PO7','O1','OZ','O2'} ; 

%vqs(:,20:80,195:350) = 0 ; 
%for i=1:size(vqs,1) ; dqs(i,:,:) = mat2gray((squeeze(vqs(i,:,:)))-imfilter((squeeze(vqs(i,:,:))),fspecial('gaussian',35,35))) ; end

th = 4 ; 
for i=1:12
   dqi = squeeze(vqs(i,:,:)) ;  
   a = dqi ; %(dqi-imfilter(dqi,fspecial('gaussian',35,25))) ; 
   za = (reshape(a,[1,551*551])) ; 
   zainds = find((za~=0)) ; zvals = (zscore(za(zainds))) ; 
   bads = find(zvals>th | zvals<-th) ; 
   za(zainds(bads)) = 0; 
   dqs(i,:,:) = reshape(za,[551,551]) ; 
end

for i=1:12 ; 
[xind,yind] = ind2sub(size(squeeze(dqs(i,:,:))),find(squeeze(dqs(i,:,:))~=0)) ; 
dqsi = squeeze(dqs(i,:,:)) ; 
vals = dqsi((squeeze(dqs(i,:,:))~=0)) ; 
[xg1,yg1] = meshgrid(1:551,1:551) ; 
itp = griddata(xind,yind,vals,xg1,yg1) ; 
idqs(i,:,:) = (itp)' ; 
end

im1 = (uint8(mat2gray(squeeze(mean(idqs(5:6,:,:),1)))*255)) ;  
im2 = (uint8(mat2gray(squeeze(mean(idqs(3:4,:,:),1)))*255)) ;  
im3 = (uint8(mat2gray(squeeze(mean(idqs(1:2,:,:),1)))*255)) ;  
expon = 2 ; 
rgbs(:,:,3) = uint8(mat2gray(im1).^expon*255) ; rgbs(:,:,2) = uint8(mat2gray(im2).^expon*255) ; rgbs(:,:,1) = uint8(mat2gray(im3).^expon*255) ;
fhandle = figure('Position',[10,-10,1000,1000]) ; 
imagesc(rgbs) ; set(gca,'XTick',[],'YTick',[])  ;
export_fig('raw_rgb.eps')
n = 65 ; coordinates = zeros(n,2) ; hold on ; 
for i=1:n
title(elecorder{i}) ; 
[x, y] = ginput(1);
coordinates(i,:) = [x, y];
plot(coordinates(:,1), coordinates(:,2),'.','Color',[1,1,0],'LineWidth',2);
text(x,y,elecorder{i},'color','w','Fontsize',12) ; 
end
hold off
export_fig('label_rgb.eps')

save_nii(make_nii(idqs),'dqs.nii.gz') ; 
pancakecoords = coordinates ; save('pancakecoords','pancakecoords') ; 
save('rgbs','rgbs') ; 
roundcoords = round(coordinates) ; 
for i=1:size(roundcoords,1) ; xgvals(i) = xg(roundcoords(i,2),roundcoords(i,1)) ; ygvals(i) = yg(roundcoords(i,2),roundcoords(i,1)) ; end
[inv_theta,inv_phi] = cart2pol(xgvals,ygvals) ; inv_phi = 1.5708-inv_phi ; 

for elec=1:size(roundcoords,1)
oz_theta = inv_theta(elec) ; oz_phi = inv_phi(elec) ; 
sumdiffs = sqrt((ftheta-oz_theta).^2+(fphi-oz_phi).^2) ;
minsumdiff_index = find(sumdiffs == min(sumdiffs),1) ; % index of the closest theta,phi point
min_rho = frho(minsumdiff_index) ; 
min_phi = fphi(minsumdiff_index) ; 
min_theta = ftheta(minsumdiff_index) ; 
[x,y,z] = sph2cart(min_theta,min_phi,min_rho) ; 
ex(elec) = x + cx ; ey(elec) = y + cy ; ez(elec) = z + cz ; 
end
elecimg = zeros(size(gsimg)) ; 
colorelecs = zeros(size(gsimg)) ; 
ex = round(ex) ; ey = round(ey) ; ez = round(ez) ; 
for coord=1:size(coordinates,1) ; 
    elecimg(ex(coord),ey(coord),ez(coord)) = 1000 ; 
    colorelecs(ex(coord),ey(coord),ez(coord)) = coord ; 
end
rute.img = imdilate(elecimg,strel(ones(3,3,3))) ; save_untouch_nii(rute,'locs.nii.gz') ; 
colorelecs = imdilate(colorelecs,strel(ones(3,3,3))) ; rute.img = colorelecs ; save_untouch_nii(rute,'color_locs.nii.gz') ; 
mricoords = [ex;ey;ez] ; save('mricoords','mricoords') ; 
save('elecorder','elecorder') ; 
rute.img = layers ; save_untouch_nii(rute,'layers.nii.gz') ; rute.img = meddilate ; save_untouch_nii(rute,'meddilate.nii.gz') ; 


