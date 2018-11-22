%%% process and hand label UTE images.
clear all ; close all ; 
cd c:/shared/lastute/ ;
mongs=dir('*') ; mongs(1:2) = [] ; 

for m=1:length(mongs) ; 
    
cd(['C:\shared\lastute\',mongs(m).name]) ; ls ; 
disp('loading raw data...') ; 
rute = load_untouch_nii('res_ute.nii.gz') ; ruteorig = double(rute.img) ; %ruteorig = ruteorig - imfilter(ruteorig,fspecial('gaussian',80,40)); 
z2 = load_untouch_nii('fnirt/ute_mask.nii.gz') ;  z2.img = medfilt3(z2.img) ; %imdilate(medfilt3(imerode(double(z2.img),strel(ones(5,5,5)))),strel(ones(5,5,5))) ; 
%outermask = load_untouch_nii('ants_ute/warpmask_outer.nii.gz') ; dilouter = imdilate(outermask.img,strel(ones(13,13,13))) ; 
%z2 = load_untouch_nii('c:/shared/ATLASES/mni_mask2.nii.gz') ; z2.img = medfilt3(double(z2.img)) ; 
%outermask = load_untouch_nii('c:/shared/ATLASES/mnimask_outer.nii.gz') ; dilouter = double(imdilate(outermask.img,strel(ones(13,13,13)))) ; 
%z2.img = medfilt3(z2.img) ; 
save_untouch_nii(z2,'finalmask.nii.gz') ; 
layer2 = imdilate(z2.img==2,strel(ones(5,5,5))) ; 

prevdil = double(z2.img==1) ; outim = double(z2.img==0) ; 

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

%rute.img = double((layers>0) .* ruteorig) ; 
%rute.img(find(z2.img==2)) = 0 ; 

dilbottom = imdilate(z2.img==2,strel(ones(12,12,12))) ; 
layers(dilbottom==1) = 0 ; 
%layers = layers.*dilouter ; 
inds = find(layers>0) ; 
%layers(inds(zscore(ruteorig(inds))<-5)) = 0 ; 
%layers(inds(zscore(ruteorig(inds))>5)) = 0 ; 

rute.img = single(layers) ; save_untouch_nii(rute,'layers.nii.gz') ; 
disp('performing pancake projection...') ; 
gsimg = double(ruteorig-imfilter(rute.img,fspecial('gaussian',61,61))) ;
% hard-coded grid variables
max_xp = 2.5 ; min_xp = -2.5 ; max_yp = 2.5 ; min_yp = -2.5 ; 
layer = layers ; 
for layer_index = 1:6
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

elecorder = {'FP1','FPZ','FP2','AF8','AF4','GND','AF3','AF7','F7','F5','F3','F1','FZ','F2','F4','F6','F8','FT10','FT8','FC6','FC4','FC2','REF','FC1','FC3','FC5','FT7','FT9',...
    'T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP10','TP8','CP6','CP4','CP2','CPZ','CP1','CP3','CP5','TP7','TP9','P7','P5','P3','P1','PZ','P2','P4','P6','P8',...
    'PO8','PO4','POZ','PO3','PO7','O1','OZ','O2'} ; 

allvqs(m,:,:,:) = vqs ; 

for i=1:size(vqs,1) ; dqs(i,:,:) = mat2gray((squeeze(vqs(i,:,:))) - imfilter((squeeze(vqs(i,:,:))),fspecial('gaussian',60,30))) ; end
th = 5 ; 
for i=1:6
   dqi = squeeze(dqs(i,:,:)) ;  
   a = dqi ; %(dqi-imfilter(dqi,fspecial('gaussian',35,25))) ; 
   za = (reshape(a,[1,551*551])) ; 
   zainds = find((za~=0)) ; zvals = (zscore(za(zainds))) ; 
   bads = find(zvals>th | zvals<-th) ; 
   za(zainds(bads)) = 0; 
   dqs(i,:,:) = reshape(za,[551,551]) ; 
end

for i=1:6 ; 
[xind,yind] = ind2sub(size(squeeze(dqs(i,:,:))),find(squeeze(dqs(i,:,:))~=0)) ; 
dqsi = squeeze(dqs(i,:,:)) ; 
vals = dqsi((squeeze(dqs(i,:,:))~=0)) ; 
[xg1,yg1] = meshgrid(1:551,1:551) ; 
itp = griddata(xind,yind,vals,xg1,yg1) ; 
idqs(i,:,:) = (itp)' ; 
end

allidqs(m,:,:,:) = idqs ; 

figure,
sf = detectSURFFeatures(squeeze(mean(idqs(3,:,:),1)),'MetricThreshold',100) ; 
imshow(squeeze(idqs(3,:,:))) ; hold on ; 
plot(selectStrongest(sf,150)) ; 

end

foregroundDetector = vision.ForegroundDetector('NumGaussians', 3, ...
    'NumTrainingFrames', 50);


vq1 = uint8(mat2gray(squeeze(allidqs(1,4,:,:)))*255) ;    
sf = detectSURFFeatures(vq1,'MetricThreshold',200,'NumOctaves',15,'NumScaleLevels',5) ; 
pts = selectStrongest(sf,200) ; 
locs = pts.Location ; 
[features1, points1] = extractFeatures(vq1, locs);

vq2 = uint8(mat2gray(squeeze(allidqs(1,4,:,:)))*255) ;    
sf = detectSURFFeatures(vq2,'MetricThreshold',200,'NumOctaves',15,'NumScaleLevels',5) ; 
pts = selectStrongest(sf,200) ; 
locs = pts.Location ; 
[features2, points2] = extractFeatures(vq2, locs);

pairs = matchFeatures(features1,features2) ; 

matched1 = points1(pairs(:,1),:) ; 
matched2 = points2(pairs(:,2),:) ; 

for i=1:14 ; 
    subplot(3,5,i) ;
    a = (squeeze(max(allvqs(i,3:end,:,:),[],2))) ; 
    imagesc(a-imfilter(a,fspecial('gaussian',25,15))) ; 
end

for i=1:14
    figure,
a = (squeeze(max(allvqs(i,3:end,:,:),[],2))) ; 
vq1 = uint8(mat2gray(a-imfilter(a,fspecial('gaussian',25,15)))*255) ; 
sf = detectSURFFeatures(vq1,'MetricThreshold',200,'NumOctaves',15,'NumScaleLevels',5) ; 
imshow(vq1) ; hold on ; 
pts = selectStrongest(sf,200) ; 
locs = pts.Location ; 
plot(locs(:,1),locs(:,2),'r.') ; 
end


