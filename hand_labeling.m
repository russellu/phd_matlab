%%% process and hand label UTE images.
%%% requires a head mask (ute_mask.nii.gz) and a 1mm isotropic UTE image
%%% requires the image processing toolbox (imdilate, imerode, imfilter, as
%%% well as some other small functions which should be included in the
%%% package, also requires the matlab nifti package by Jimmy Shen
%%% outputs a bunch of stuff, but the original hand-labeled coordinates are
%%% saved in mricoords_1, and are in "voxel space" of the UTE which was
%%% input to the program, ie each coordinate is a voxel offset along that
%%% axis (x,y or z). 
% author Russell Butler russell.buttler@usherbrooke.ca

clear all ; close all ; 

% paths to relevant data and output 
path_ute = 'c:/shared/lastute/russ/res_ute.nii.gz' ; 
path_mask = 'c:/shared/lastute/russ/fnirt/ute_mask.nii.gz' ; 
path_output = 'c:/shared/ute_output/russ' ;
mkdir(path_output) ; cd(path_output) ; 

% the order in which to label the electrodes. after each mouse click, the
% program will move to the next electrode, ie, each mouse click labels an
% electrode in the sequence "elecorder" (wait for the crosshair)
% it is recommended that you put a gel capsule on the cap, before the
% experiment so you can be 100% sure which side is left and which is right,
% as this might get swapped depending on the converter you use
elecorder = {'FP1','FPZ','FP2','AF8','AF4','GND','AF3','AF7','F7','F5','F3','F1','FZ','F2','F4','F6','F8','FT10','FT8','FC6','FC4','FC2','REF','FC1','FC3','FC5','FT7','FT9',...
    'T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP10','TP8','CP6','CP4','CP2','CPZ','CP1','CP3','CP5','TP7','TP9','P7','P5','P3','P1','PZ','P2','P4','P6','P8',...
    'PO8','PO4','POZ','PO3','PO7','O1','OZ','O2'} ; 

% get the raw UTE, and the mask
disp('loading raw data...') ; 
rute = load_untouch_nii(path_ute) ; ruteorig = double(rute.img) ; % raw ute
mask = load_untouch_nii(path_mask) ;  % mask

% save UTE and mask in the output directory, for segmentation script
save_untouch_nii(rute,'res_ute.nii.gz') ; 
save_untouch_nii(mask,'finalmask.nii.gz') ; 

% save the final mask in the output directory
save_untouch_nii(mask,'finalmask.nii.gz') ; 

% create the layers surrounding the head mask (12 layers)
layer2 = imdilate(mask.img==2,strel(ones(5,5,5))) ; 
prevdil = double(mask.img==1) ; outim = double(mask.img==0) ; 
disp('performing iterative dilation...') ; 
intensitylayers = zeros(size(prevdil)) ;
for i=1:12 ; 
    dilmaski = imdilate(prevdil,strel(ones(3,3,3))) .* outim ; 
    clayer = (dilmaski - prevdil) > 0 ; 
    intensitylayers(clayer==1) = i ; 
    if i==1
        prevdil = dilmaski ;
    else
        prevdil = (prevdil + dilmaski) > 0 ; 
    end
    intensitylayers = intensitylayers.*(~layer2) ; 
end
disp('computing final layers...') ; 
layers = intensitylayers ;
dilbottom = imdilate(mask.img==2,strel(ones(12,12,12))) ; 
layers(dilbottom==1) = 0 ; 
inds = find(layers>0) ; 

% save the layers
rute.img = single(layers) ; save_untouch_nii(rute,'layers.nii.gz') ; 

% do the pancake projection
disp('performing pancake projection...') ; 
gsimg = double(ruteorig-imfilter(rute.img,fspecial('gaussian',61,61))) ; % first, take the gradient of the raw UTE (so the electrodes show up more clearly)
max_xp = 2.5 ; min_xp = -2.5 ; max_yp = 2.5 ; min_yp = -2.5 ; % hard-coded grid variables
for layer_index = 1:6
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

for i=1:size(vqs,1) ; 
    dqs(i,:,:) = mat2gray((squeeze(vqs(i,:,:))) - imfilter((squeeze(vqs(i,:,:))),fspecial('gaussian',60,30))) ;
end

% interpolate super high intensity pixels (z>5), for improved image contrast
% first find the bad pixels
th = 5 ; 
for i=1:6
   dqi = squeeze(dqs(i,:,:)) ;  
   za = (reshape(dqi,[1,551*551])) ; 
   zainds = find((za~=0)) ; zvals = zscore(za(zainds)) ; 
   bads = find(zvals>th | zvals<-th) ; 
   za(zainds(bads)) = 0; 
   dqs(i,:,:) = reshape(za,[551,551]) ; 
end
% then do the interpolation
for i=1:6 ; 
    [xind,yind] = ind2sub(size(squeeze(dqs(i,:,:))),find(squeeze(dqs(i,:,:))~=0)) ; 
    dqsi = squeeze(dqs(i,:,:)) ; 
    vals = dqsi((squeeze(dqs(i,:,:))~=0)) ; 
    [xg1,yg1] = meshgrid(1:551,1:551) ; 
    itp = griddata(xind,yind,vals,xg1,yg1) ; 
    idqs(i,:,:) = (itp)' ; 
end

% create the RGB pancake
im1 = (uint8(mat2gray(squeeze(mean(idqs(5:6,:,:),1)))*255)) ;  
im2 = (uint8(mat2gray(squeeze(mean(idqs(3:4,:,:),1)))*255)) ;  
im3 = (uint8(mat2gray(squeeze(mean(idqs(1:2,:,:),1)))*255)) ;  
expon = 1.5 ; % to accentuate higher value pixels
rgbs(:,:,3) = uint8(mat2gray(im1).^expon*255) ; rgbs(:,:,2) = uint8(mat2gray(im2).^expon*255) ; rgbs(:,:,1) = uint8(mat2gray(im3).^expon*255) ; 
imwrite(rgbs,'rgbs.png') ;
save_nii(make_nii(idqs(1:6,:,:)),'idqs.nii.gz') ; 

maskvq = zeros(size(vqs)) ;
maskvq(vqs~=0) = 1 ;

save_nii(make_nii(maskvq(1:5,:,:)),'maskvq.nii.gz') ; 

fhandle = figure('Position',[10,-10,1000,1000]) ; 
imagesc(rgbs) ; 
n = 65 ; coordinates = zeros(n,2) ; hold on ; 
for i=1:n
title(elecorder{i}) ; 
[x, y] = ginput(1);
coordinates(i,:) = [x, y];
plot(coordinates(:,1), coordinates(:,2),'.','Color',[1,1,0],'LineWidth',2);
text(x,y,elecorder{i},'color','w','Fontsize',12) ; 
end
hold off

% randomize the coordinates, and invert them to 3d MRI space (the original
% hand-labeled coordinates are in mricoords_1)
for mricoordn=1:30 ; 
    if mricoordn==1 ; % save the original coordinates
        randcoords(:,1) = coordinates(:,1) ; 
        randcoords(:,2) = coordinates(:,2) ; 
    else
        randcoords(:,1) = coordinates(:,1) + (rand(1,65)'-.5)*10 ; % random offset (produces ~1.5mm standard deviation in 3d MRI space)
        randcoords(:,2) = coordinates(:,2) + (rand(1,65)'-.5)*10 ; 
    end
    roundcoords = round(randcoords) ; 
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
    mricoords = [ex;ey;ez] ; save(['mricoords_',num2str(mricoordn)],'mricoords') ; 
end

