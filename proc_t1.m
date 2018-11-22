clear all ; close all ; 
cd c:/shared/lastute/ ;
mongs=dir('*') ; mongs(1:2) = [] ; 
mricoordn = 4; 

for m=1:length(mongs) ; 
cd(['c:/shared/lastute/',mongs(m).name]) ; ls ; 
t1 = load_untouch_nii('t1_in_ute.nii.gz') ; 
mask = load_untouch_nii('finalmask.nii.gz') ; mask = medfilt3(mask.img) ; 

layer = mask - imerode(mask,strel(ones(3,3,3))) ; 
dilbottom = imdilate(mask==2,strel(ones(12,12,12))) ; 
layer(dilbottom==1) = 0 ; 

% hard-coded grid variables
max_xp = 2.5 ; min_xp = -2.5 ; max_yp = 2.5 ; min_yp = -2.5 ; 
gsimg = double(t1.img) ; layers = layer ; 
clear vqs 
for layer_index = 1
    layerinds = find(layers==layer_index) ; 
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

vqs = squeeze(uint8(mat2gray(vqs)*255)) ; vqs(vqs==0) = 255 ; 

elecorder = {'FP1','FPZ','FP2','AF8','AF4','GND','AF3','AF7','F7','F5','F3','F1','FZ','F2','F4','F6','F8','FT10','FT8','FC6','FC4','FC2','REF','FC1','FC3','FC5','FT7','FT9',...
    'T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP10','TP8','CP6','CP4','CP2','CPZ','CP1','CP3','CP5','TP7','TP9','P7','P5','P3','P1','PZ','P2','P4','P6','P8',...
    'PO8','PO4','POZ','PO3','PO7','O1','OZ','O2'} ; 

figure,imshow(vqs) ; 

allvqs(m,:,:) = vqs ;
%{
fhandle = figure('Position',[10,-10,1000,1000]) ; 
imshow(vqs) ; %set(gca,'XTick',[],'YTick',[])  ;
n = 65 ; coordinates = zeros(n,2) ; hold on ; 
for i=1:n
title(elecorder{i}) ; 
[x, y] = ginput(1);
coordinates(i,:) = [x, y];
plot(coordinates(:,1), coordinates(:,2),'.','Color',[1,1,0],'LineWidth',2);
text(x,y,elecorder{i},'color','w','Fontsize',12) ; 
end
hold off
%export_fig('label_rgb.eps')

    
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
colorelecs = zeros(size(gsimg)) ; 
ex = round(ex) ; ey = round(ey) ; ez = round(ez) ; 
for coord=1:size(coordinates,1) ; 
    elecimg(ex(coord),ey(coord),ez(coord)) = 1000 ; 
    colorelecs(ex(coord),ey(coord),ez(coord)) = coord ; 
end
mricoords = [ex;ey;ez] ; save(['gel_mricoords_',num2str(mricoordn)],'mricoords') ; 
%}
end



for i=1:14 ; subplottight(4,4,i) ; imagesc(squeeze(allvqs(i,:,:))) ; colormap gray ; text(20,20,['S ',num2str(i)]) ; set(gca,'XTick',[],'YTick',[]) ; end



for i=1:14 ; subplottight(4,4,i) ; imagesc(squeeze(allrgbs(i,:,:,:))) ; text(20,20,['S ',num2str(i)]) ; set(gca,'XTick',[],'YTick',[]) ; end








