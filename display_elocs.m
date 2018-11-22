%%%%% converting the pancake properly to and from a volume
clear all ; close all ; 
cd c:/shared/ute/ ; 

subs = dir('*') ; subs(1:2) = [] ; 
for sb=1:length(subs) ; 
    cd(['c:/shared/ute/',subs(sb).name]) ; 
    intlayers = load_untouch_nii('intensity_layers.nii.gz') ; 
    gs = load_untouch_nii('res_ute.nii.gz') ; 
    gsimg = double(gs.img-imfilter(gs.img,fspecial('gaussian',41,41))) ;
    % hard-coded grid variables
    max_xp = 2.5 ; min_xp = -2.5 ; max_yp = 2.5 ; min_yp = -2.5 ; 
    layer = intlayers.img>0 ;  
    for layer_index = 1:8
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

    %vqs(:,135:160,150:180) = 0 ; 
    for i=1:size(vqs,1) ; dqs(i,:,:) = (squeeze(vqs(i,:,:)))-imfilter((squeeze(vqs(i,:,:))),fspecial('gaussian',40,40)) ; end
    rgbs(:,:,3) = uint8(mat2gray(squeeze(mean(dqs(4,:,:),1)))*255) ;  
    rgbs(:,:,2) = uint8(mat2gray(squeeze(mean(dqs(3,:,:),1)))*255) ;  
    rgbs(:,:,1) = uint8(mat2gray(squeeze(mean(dqs(1:2,:,:),1)))*255) ;  

    %fhandle = figure('Position',[10,-10,1000,1000]) ; 
    subplot(1,8,sb) ; 
    imagesc(rgbs) ; set(gca,'XTickLabel',[],'YTickLabel',[]) ; 
end