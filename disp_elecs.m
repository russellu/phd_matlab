%%% process and hand label UTE images.
clear all ; close all ;
subs = {'alex','biz','cloud','dave','felix','jeremie','karl','nic','pierre','russ','sukh','terry'} ; 
for subj=6%:length(subs) ; 
name = subs{subj} ; 
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

rute.img = layers ; save_untouch_nii(rute,'layers.nii.gz') ; 

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
%fhandle = figure('Position',[10,-10,1000,1000]) ; 
subplot(3,4,subj) ; 
imagesc(rgbs) ; set(gca,'XTick',[],'YTick',[])  ;

allrgbs(subj,:,:,:) = rgbs ; 
end


for i=1:12 ; 
    h=tight_subplot(3,4,.01,[.1 .01],[.01 .01]) ;
   
end
for i=1:12 ; axes(h(i)) ; imshow(squeeze(allrgbs(i,:,:,:))) ; end

cd c:/shared/utef/jeremie ; ls
layers = load_untouch_nii('layers.nii.gz') ; 

layers(layers>6) = 0 ; 
rgblayer = zeros(240,240,170,3) ; rgblayer(:,:,:,1) = layers>0 & layers<3 ; rgblayer(:,:,:,2) = layers>2 & layers<5 ; rgblayer(:,:,:,3) = layers>4 & layers<7 ; 
plotoverlayIntensity2D(squeeze(ruteorig(:,:,75)),(mat2gray(squeeze(layers(:,:,75)))>0)*.5,(squeeze(layers(:,:,75))>0),270) ;

anat = squeeze(ruteorig(:,:,75)) ; rotangle = 270 ; 
heatmap = squeeze(layers(:,:,75)>0) ; 
img = imrotate(uint8(mat2gray(anat)*255),rotangle) ; img(isnan(img)) = 0 ; img(isinf(img)) = 0 ; 
im2 = imrotate(uint8(mat2gray(heatmap)*255),rotangle) ; im2(isnan(im2)) = 0 ; im2(isinf(im2)) = 0 ;
I = imrotate(uint8(mat2gray(layers(:,:,75))*255),270) ;  
rgb = ind2rgb(gray2ind(I,255),hot(255)) ;
%for i=1:3
%    rlayers(:,:,i) = imrotate(squeeze(rgblayer(:,:,75,i)),rotangle) ; 
%end
%rgb = squeeze(rlayers) ; 
imshow(rgb) ; 
imshow(img, 'InitialMag', 'fit') ;
hold on ; h = imshow(rgb) ; hold off ; 
set(h, 'AlphaData', imrotate((squeeze(heatmap)),rotangle)) ;

