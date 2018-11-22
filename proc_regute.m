%%% process and hand label UTE images.


clear all ; close all ; 
cd c:/shared/lastute/ ;
mongs=dir('*') ; mongs(1:2) = [] ; 

for m=1%:length(mongs) ; 
    
%cd(['C:\shared\lastute\',mongs(m).name]) ; ls ; 
cd(['C:\shared\lastute\genevieve']) ; ls ; 

disp('loading raw data...') ; 
rute = load_untouch_nii('res_ute.nii.gz') ; ruteorig = double(rute.img) ; %ruteorig = ruteorig - imfilter(ruteorig,fspecial('gaussian',80,40)); 
z2 = load_untouch_nii('mask.nii.gz') ;  z2.img = medfilt3(z2.img) ; %imdilate(medfilt3(imerode(double(z2.img),strel(ones(5,5,5)))),strel(ones(5,5,5))) ; 
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

im1 = (uint8(mat2gray(squeeze(mean(idqs(5,:,:),1)))*255)) ;  
im2 = (uint8(mat2gray(squeeze(mean(idqs(3:4,:,:),1)))*255)) ;  
im3 = (uint8(mat2gray(squeeze(mean(idqs(1:2,:,:),1)))*255)) ;  
expon = 1.5 ; 
rgbs(:,:,3) = uint8(mat2gray(im1).^expon*255) ; rgbs(:,:,2) = uint8(mat2gray(im2).^expon*255) ; rgbs(:,:,1) = uint8(mat2gray(im3).^expon*255) ; 
%figure,
%imagesc(rgbs) ; %set(gca,'XTick',[],'YTick',[])  ;
imwrite(rgbs,'rgbs.png') ;
allrgbs(m,:,:,:) = rgbs ; 
save_nii(make_nii(idqs(1:5,:,:)),'idqs.nii.gz') ; 

%
maskvq = zeros(size(vqs)) ;
maskvq(vqs~=0) = 1 ;

save_nii(make_nii(maskvq(1:5,:,:)),'maskvq.nii.gz') ; 

fhandle = figure('Position',[10,-10,1000,1000]) ; 
imagesc(rgbs) ; %set(gca,'XTick',[],'YTick',[])  ;
n = 65 ; coordinates = zeros(n,2) ; hold on ; 
for i=1:n
title(elecorder{i}) ; 
[x, y] = ginput(1);
coordinates(i,:) = [x, y];
plot(coordinates(:,1), coordinates(:,2),'.','Color',[1,1,0],'LineWidth',2);
text(x,y,elecorder{i},'color','w','Fontsize',12) ; 
end
hold off

for mricoordn=1:30 ; 
    if mricoordn==1;
        randcoords(:,1) = coordinates(:,1) ; 
        randcoords(:,2) = coordinates(:,2) ; 
    else
        randcoords(:,1) = coordinates(:,1) + (rand(1,65)'-.5)*10 ; 
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

end


%{
fhandle = figure('Position',[10,-10,1000,1000]) ; 
imagesc(rgbs) ; %set(gca,'XTick',[],'YTick',[])  ;
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

pancakecoords = coordinates ; save('pancakecoords','pancakecoords') ; 
save('rgbs','rgbs') ; 
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
rute.img = imdilate(elecimg,strel(ones(3,3,3))) ; save_untouch_nii(rute,'locs.nii.gz') ; 
colorelecs = imdilate(colorelecs,strel(ones(3,3,3))) ; rute.img = colorelecs ; save_untouch_nii(rute,'color_locs.nii.gz') ; 
mricoords = [ex;ey;ez] ; save('mricoords','mricoords') ; 
save('elecorder','elecorder') ; 
rute.img = layers ; save_untouch_nii(rute,'layers.nii.gz') ; %rute.img = meddilate ; save_untouch_nii(rute,'meddilate.nii.gz') ; 

%}

%{
ref = load_untouch_nii('res_ute.nii.gz') ; 
cd c:/shared/lastute/biz ; ls 
mcoords = dir('mricoords_*') ; 
for m=1:length(mcoords)
    coordi = load(mcoords(m).name) ; coordi = coordi.mricoords ; 
    refimg = zeros(size(ref.img)) ; 
    for j=1:size(coordi,2)
        refimg(coordi(1,j),coordi(2,j),coordi(3,j)) = 1 ; 
    end
    ref.img = imdilate(refimg,strel(ones(3,3,3))) ; 
    save_untouch_nii(ref,['rndcoords_',num2str(m),'.nii.gz']) ; 
end



for i=1:14
    subplottight(4,4,i) ;
    imagesc(squeeze(allrgbs(i,:,:,:))) ;
    set(gca,'XTick',[],'YTick',[]) ; 
    text(20,20,['S',num2str(i)],'Color',[1,1,1]) ; 
end

%}



