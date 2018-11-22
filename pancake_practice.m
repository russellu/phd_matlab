%%%%% converting the pancake properly to and from a volume
clear all ; close all ; 
cd c:/shared/tests2/ ; ls 
intlayers = load_untouch_nii('intensity_layers.nii.gz') ; 

gs = load_untouch_nii('res_3.nii.gz') ; 
gsimg = double(gs.img-imfilter(gs.img,fspecial('gaussian',41,41))) ;
% hard-coded grid variables
max_xp = 2.5 ; min_xp = -2.5 ; max_yp = 2.5 ; min_yp = -2.5 ; 
layer = intlayers.img>0 ;  
for layer_index = 1:12
    layerinds = find(intlayers.img==layer_index) ; 
    [lx,ly,lz] = ind2sub(size(layer),layerinds) ; 
    [cx,cy,cz] = centmass3(layer) ; 
    xlayer_diffs = lx-cx ; ylayer_diffs = ly-cy ; zlayer_diffs = lz-cz ; 
    [theta,phi,rho] = cart2sph(xlayer_diffs,ylayer_diffs,zlayer_diffs) ; 
    if layer_index==1 ; ftheta = theta ; fphi = phi ; frho = rho ; end
    theta_brain = zeros(size(layer)) ; phi_brain = zeros(size(layer)) ; rho_brain = zeros(size(layer)) ;
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
for i=1:size(vqs,1) ; dqs(i,:,:) = mat2gray((squeeze(vqs(i,:,:)))-imfilter((squeeze(vqs(i,:,:))),fspecial('gaussian',40,40))) ; end
rgbs(:,:,3) = uint8(mat2gray(squeeze(mean(dqs(7:8,:,:),1)))*255) ;  
rgbs(:,:,2) = uint8(mat2gray(squeeze(mean(dqs(5:6,:,:),1)))*255) ;  
rgbs(:,:,1) = uint8(mat2gray(squeeze(mean(dqs(3:4,:,:),1)))*255) ;  

fhandle = figure('Position',[10,-10,1000,1000]) ; 
imagesc(rgbs) ; 
n = 65 ; coordinates = zeros(n,2) ; hold on ; 
for i=1:n
title(elecorder{i}) ; 
[x, y] = ginput(1);
coordinates(i,:) = [x, y];
plot(coordinates(:,1), coordinates(:,2),'+','Color',[1,1,0],'LineWidth',2);
end
hold off

pancakecoords = coordinates ; save('pancakecoords','pancakecoords') ; 

roundcoords = round(coordinates) ; 
for i=1:size(roundcoords,1) ; xgvals(i) = xg(roundcoords(i,2),roundcoords(i,1)) ; ygvals(i) = yg(roundcoords(i,2),roundcoords(i,1)) ; end
[inv_theta,inv_phi] = cart2pol(xgvals,ygvals) ; inv_phi = 1.5708-inv_phi ; 
%%% now that you have inverted the coordinates, you need to find the
%%% closest 3D point that matches these
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
ex = round(ex) ; ey = round(ey) ; ez = round(ez) ; 
for coord=1:size(coordinates,1) ; 
    elecimg(ex(coord),ey(coord),ez(coord)) = 1000 ; 
end
gs.img = elecimg ; save_untouch_nii(gs,'locs.nii.gz') ; 
mricoords = [ex;ey;ez] ; save('mricoords','mricoords') ; 

%{
cd C:\shared\elec_stds ; 
ls 
hres = load('hr_carts.mat') ; hres = hres.carts ; stdh = squeeze(std(hres,0,2)) ; 
mres = load('mr_carts.mat') ; mres = mres.carts ; stdm = squeeze(std(mres,0,2)) ; 
lres = load('lr_carts.mat') ; lres = lres.carts ; stdl = squeeze(std(lres,0,2)) ; 
bchart = [mean(stdl,1); mean(stdm,1); mean(stdh,1)] ; 
bar(bchart') ; 
set(gca,'XTick',1:65,'XTickLabel',elecorder) ; 
rotateticklabel(gca,45) ; ylabel('std(mm)') ; hline(1,'k') ; 
reslabs = {'2.5mm isotropic','1.875x1.875x2mm','1.5mm isotropic'} ; 
barwitherr(squeeze(std(bchart,0,2)),squeeze(mean(bchart,2))) ;
set(gca,'XTickLabel',reslabs) ; ylabel('std(mm)') ; 
allcoords(5,:,:) = coordinates ; 
colors = {[1,0,0],[1,1,0],[0,1,0],[0,1,1],[1,1,1]} ; 
fhandle = figure('Position',[10,-10,1000,1000]) ; 
imagesc(rgbs) ; hold on ; 
for i=1:size(allcoords,1) ; plot(allcoords(i,:,1),allcoords(i,:,2),'+','Color',colors{i},'LineWidth',2) ; end 
clear ex ey ez allex alley allez
for ac = 1:size(allcoords,1) ; 
coordinates = squeeze(allcoords(ac,:,:)) ; 
roundcoords = round(coordinates) ; 
for i=1:size(roundcoords,1) ; xgvals(i) = xg(roundcoords(i,2),roundcoords(i,1)) ; ygvals(i) = yg(roundcoords(i,2),roundcoords(i,1)) ; end
[inv_theta,inv_phi] = cart2pol(xgvals,ygvals) ; inv_phi = 1.5708-inv_phi ; 
%%% now that you have inverted the coordinates, you need to find the
%%% closest 3D point that matches these
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
allex(ac,:) = (ex) ; alley(ac,:) = (ey) ; allez(ac,:) = (ez) ; 
%for coord=1:size(coordinates,1) ; 
%    elecimg(ex(coord)-1:ex(coord)+1,ey(coord)-1:ey(coord)+1,ez(coord)-1:ez(coord)+1) = 1000 ; 
%end
%gs.img = elecimg ; save_untouch_nii(gs,['test_',num2str(ac),'.nii.gz']) ; 
end
%mean_std = std(allez,0,1) + std(allex,0,1) + std(alley,0,1) / 3 ;
%carts(1,:,:) = allex ; carts(2,:,:) = alley ; carts(3,:,:) = allez ;  
%cd('c:/shared/elec_stds') ; save('mr_carts','carts') ; 
%}