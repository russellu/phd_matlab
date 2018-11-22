%%% process and hand label UTE images.
clear all ; close all ; 
cd('C:\shared\regute\') ; ls ; 
disp('loading raw data...') ; 
rute = load_untouch_nii('dof12_felix.nii.gz') ; ruteorig = rute.img ; 
rutime = rute.img ; 
%{
[~,~,zg] = ndgrid(1:size(rutime,1),1:size(rutime,2),1:size(rutime,3)) ; 
rutime = rutime.*(zg.^1.5) ;% bias correction (z)
%rutime = (rutime-imfilter(rutime,fspecial('gaussian',15,15))) ;
%gradu = rutime - imfilter(rutime,fspecial('gaussian',35,35)) ;
disp('performing bilateral intensity normalization') ; 
sym = symmetrize(rute) ;
%sym = rutime ; 

resute = sym(1:2:end,1:2:end,1:2:end) ; 
disp('segmenting using k-means...') ; 
[kc,km] = kmeans(uint8(mat2gray(reshape(resute,[1,numel(resute)]))*255),2) ; 
kimg = reshape(km,size(resute)) > 1 ; 
padamt = 60 ; 
disp('padding image...') ; 
padded = pad3d(kimg,padamt) ; 
padk = pad3d(kimg,padamt) ; 
disp('calculating center of mass...') ; 
[cx,cy,cz] = centmass3(abs(padded)) ; [xg,yg,zg] = ndgrid(-cx:size(padded,1)-cx-1,-cy:size(padded,2)-cy-1,-cz:size(padded,3)-cz-1) ; 
disp('transforming to spherical coordinates...') ; 
[th,phi,rho] = cart2sph(xg,yg,zg) ; 
rhomask = rho.*(padk) ; 
disp('computing bounding sphere...') ; 
sphere = sqrt(xg.^2 + yg.^2 + zg.^2) < 80 ; rhomask = rhomask.*sphere ; 
shell = sphere-imerode(sphere,strel(ones(3,3,3))) ; shellinds = find(shell==1) ; 
[sx,sy,sz] = ind2sub(size(shell),shellinds) ;
diffx = sx-cx ; diffy = sy-cy ; diffz = sz-cz ; difflengths = sqrt(diffx.^2 + diffy.^2 + diffz.^2) ; 
udiffx = diffx./difflengths ; udiffy = diffy./difflengths ; udiffz = diffz./difflengths ;
disp('creating projection lines...') ; 
lines = repmat((1:round(mean(difflengths)))',[1,length(difflengths)]) ; 
clear xlines ylines zlines
for i=1:size(lines,1) ; 
    xlines(i,:) = floor(cx - lines(i,:).*udiffx')+1 ; 
    ylines(i,:) = floor(cy - lines(i,:).*udiffy')+1 ; 
    zlines(i,:) = floor(cz - lines(i,:).*udiffz')+1 ; 
end
disp('getting line indices...') ; 
clear inds vals
for i=1:size(lines,1)
    inds(i,:) = sub2ind(size(shell),xlines(i,:),ylines(i,:),zlines(i,:)) ; 
    vals(i,:) = rhomask(inds(i,:)) ; 
end
disp('creating head boundary...') ; 
[maxvals,maxinds] = max(vals,[],1) ; 
zbrain = zeros(size(padded)) ; 
for i=1:size(maxinds,2)
    indi = inds(maxinds(i),i) ;
    zbrain(indi) = 1 ; 
end
newbrain = zbrain(padamt:size(zbrain,1)-padamt-1,padamt:size(zbrain,2)-padamt-1,padamt:size(zbrain,3)-padamt-1) ; 
%}
disp('filling head with ray tracing...') ; 
%%%% do a region growing approach (using ray tracing)
%{
slicebrain = zeros(size(newbrain)) ; 
for i=1:size(newbrain,3) 
     slicei = squeeze(newbrain(:,:,i)) ; 
     for j=1:size(slicei,1) % for all rows
         filling = true ; fillcount = 1 ; 
         for k=1:size(slicei,2) % for all columns
             if slicei(j,k)==1 
                 fillcount = fillcount + 1 ;                 
             end
             if fillcount > 1 && slicei(j,k) == 0
                fillcount = 1 ; filling = ~filling ;  
             end
             if filling 
                 slicei(j,k) = 1 ;                  
             end                           
         end
     end   
     slicebrain(:,:,i) = slicei ; 
end
slicebrain2 = zeros(size(newbrain)) ; 
for i=1:size(newbrain,3) 
     slicei = squeeze(newbrain(:,:,i)) ; 
     for j=1:size(slicei,2) % for all rows
         filling = true ; fillcount = 1 ; 
         for k=1:size(slicei,1) % for all columns
             if slicei(k,j)==1 
                 fillcount = fillcount + 1 ;                 
             end
             if fillcount > 1 && slicei(j,k) == 0
                fillcount = 1 ; filling = ~filling ;  
             end
             if filling 
                 slicei(k,j) = 1 ;                  
             end                           
         end
     end   
     slicebrain2(:,:,i) = slicei ; 
end
slicebrain3 = zeros(size(newbrain)) ; 
for i=1:size(newbrain,3) 
     slicei = squeeze(newbrain(:,:,i)) ; 
     for j=size(slicei,2):-1:1 % for all rows
         filling = true ; fillcount = 1 ; 
         for k=size(slicei,1):-1:1 % for all columns
             if slicei(k,j)==1 
                 fillcount = fillcount + 1 ;                 
             end
             if fillcount > 1 && slicei(j,k) == 0
                fillcount = 1 ; filling = ~filling ;  
             end
             if filling 
                 slicei(k,j) = 1 ;                  
             end                           
         end
     end   
     slicebrain3(:,:,i) = slicei ; 
end
slicebrain4 = zeros(size(newbrain)) ; 
for i=1:size(newbrain,3) 
     slicei = squeeze(newbrain(:,:,i)) ; 
     for j=size(slicei,1):-1:1 % for all rows
         filling = true ; fillcount = 1 ; 
         for k=size(slicei,2):-1:1 % for all columns
             if slicei(j,k)==1 
                 fillcount = fillcount + 1 ;                 
             end
             if fillcount > 1 && slicei(j,k) == 0
                fillcount = 1 ; filling = ~filling ;  
             end
             if filling 
                 slicei(j,k) = 1 ;                  
             end                           
         end
     end   
     slicebrain4(:,:,i) = slicei ; 
end

sumslice = slicebrain+slicebrain2+slicebrain3+slicebrain4 ; 
binslice = medfilt3(sumslice)> 2 ; 
disp('segmenting head from background...') ; 
bw = bwconncomp(binslice) ; clear lengths ; 
pixlist = bw.PixelIdxList ; for i=1:length(pixlist) ; lengths(i) = length(pixlist{i}) ; end
compimg = zeros(size(slicebrain)) ; compimg(pixlist{find(lengths==max(lengths))}) = 1 ; 
compimg = ~compimg ; 
bw = bwconncomp(compimg) ; clear lengths ; 
pixlist = bw.PixelIdxList ; for i=1:length(pixlist) ; lengths(i) = length(pixlist{i}) ; end
compimg = zeros(size(slicebrain)) ; compimg(pixlist{find(lengths==max(lengths))}) = 1 ; 
disp('performing morphological operations to smooth head surface...erosion...') ; 
mederode = medfilt3(imerode(medfilt3(imerode(compimg,strel(ones(3,3,3)))),strel(ones(3,3,3)))) ; 
disp('performing morphological operations to smooth head surface...dilation...') ; 
meddilate = medfilt3(imdilate(medfilt3(imdilate(mederode,strel(ones(3,3,3)))),strel(ones(3,3,3)))) ; 
disp('head segmentation complete.') ; 

disp('upsampling dilation mask...') ; 
[xg1,yg1,zg1] = ndgrid(1:size(meddilate,1),1:size(meddilate,2),1:size(meddilate,3)) ; 
[xg2,yg2,zg2] = ndgrid(1:.5:size(meddilate,1)+.5,1:.5:size(meddilate,2)+.5,1:.5:size(meddilate,3)+.5) ; 
p = [2,1,3] ; xg1 = permute(xg1,p) ; yg1 = permute(yg1,p) ; zg1 = permute(zg1,p) ; medv = permute(meddilate,p) ; 
vq = interp3(xg1,yg1,zg1,medv,xg2,yg2,zg2,'cubic') ; 
%}
%{
meddilate = medfilt3(vq>0) ; 
resid = ((meddilate==0).*rutime) ; 
meddilate = (meddilate + (resid>mean(mean(mean(resid)))*15)) ; 
struc3 = strel(ones(5,5,5)) ; 
meddilate = imerode(imerode(imerode(imerode(meddilate,struc3),struc3),struc3),struc3) ; 
meddilate = medfilt3(meddilate) ; 
meddilate = imdilate(imdilate(imdilate(imdilate(meddilate,struc3),struc3),struc3),struc3) ; 
%}
utemask = load_untouch_nii('maskinv.nii.gz') ; utemask = utemask.img ; meddilate = double(utemask) ; 

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
layers = intensitylayers.*(meddilate==0) ; layers(:,:,1:4) = 0 ; 
%disp('creating skin surface...') ; 
%skinsurf = (imdilate(meddilate>0,strel(ones(3,3,3))) .* maskim) ;


disp('performing pancake projection...') ; 
gsimg = double(ruteorig-imfilter(rute.img,fspecial('gaussian',41,41))) ;
% hard-coded grid variables
max_xp = 2.5 ; min_xp = -2.5 ; max_yp = 2.5 ; min_yp = -2.5 ; 
badlayers = find(gsimg>(mean(mean(mean(gsimg)))*150)) ; layers(badlayers) = 0 ; layers = medfilt3(layers) ; 
border = (layers>0) .* ruteorig ;
borderdiff = (abs(border-imfilter(border,fspecial('gaussian',3,3)))) ; 
diffinds = find(borderdiff~=0) ; diffvals = borderdiff(diffinds) ; 
layers2 = layers ; baddiffs = find(zscore(diffvals)>2) ; baddiffinds = diffinds(baddiffs) ; 
layers(baddiffinds) = 0 ; 
layer = layers>0 ;  
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
for i=1:size(vqs,1) ; dqs(i,:,:) = mat2gray((squeeze(vqs(i,:,:)))-imfilter((squeeze(vqs(i,:,:))),fspecial('gaussian',40,40))) ; end
zt = 12 ; 
im1 = (uint8(mat2gray(squeeze(mean(dqs(5:8,:,:),1)))*255)) ;  
bads = find(abs(zscore(double(reshape(im1,[1,numel(im1)]))))>zt) ; im1(bads) = mean(mean(im1)) ; 
im2 = (uint8(mat2gray(squeeze(mean(dqs(3:4,:,:),1)))*255)) ;  
bads = find(abs(zscore(double(reshape(im2,[1,numel(im2)]))))>zt) ; im2(bads) = mean(mean(im2)) ; 
im3 = (uint8(mat2gray(squeeze(mean(dqs(1:2,:,:),1)))*255)) ;  
bads = find(abs(zscore(double(reshape(im3,[1,numel(im3)]))))>zt) ; im3(bads) = mean(mean(im3)) ; 
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

save_nii(make_nii(dqs),'dqs.nii.gz') ; 
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
rute.img = layers ; save_untouch_nii(rute,'layers.nii.gz') ; rute.img = meddilate ; save_untouch_nii(rute,'meddilate.nii.gz') ; 


