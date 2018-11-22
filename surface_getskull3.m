clear all ; close all ; 
cd c:/shared/UTE ; ls 
%s = load_untouch_nii('res_lin.nii.gz') ; 
s = load_untouch_nii('res_lin.nii.gz') ; 
simg = double(s.img) ; 
gs = double(simg-imfilter(simg,fspecial('gaussian',65,300))) ; 

% k means on every slice individually in 3 separate directions x,y,z
for i=1:240 ; [~,ksx(i,:,:)] = kmeans(round(mat2gray(squeeze(simg(i,:,:)))*255),2) ; end 
for i=1:240 ; [~,ksy(:,i,:)] = kmeans(round(mat2gray(squeeze(simg(:,i,:)))*255),2) ; end
for i=1:150 ; [~,ksz(:,:,i)] = kmeans(round(mat2gray(squeeze(simg(:,:,i)))*255),2) ; end

% trim the slices (doesn't help much)
fst = load_untouch_nii('fast/fast_seg_2.nii.gz') ; bw = bwconncomp(fst.img) ; remks = zeros(size(fst.img)) ; remks(bw.PixelIdxList{1}) = 1 ; 
%{
ks = ((ksy==2) .* (ksx==2) .* (ksz==2)) ;
remks = zeros(size(ks)) ;
for i=1:size(ks,3) ; 
    bw = bwconncomp(squeeze(ks(:,:,i))) ; plist = bw.PixelIdxList ; 
    maxsize = 1 ; maxind = 1 ; 
    for p=1:length(plist) ; if length(plist{p}) > maxsize ; maxind = p ; maxsize = length(plist{p}) ; end ; end 
    slice = zeros(size(ks,1),size(ks,2)) ; 
    if ~isempty(plist) ; slice(bw.PixelIdxList{maxind}) = 1 ; end
    remks(:,:,i) = slice ; 
end 
for i=1:size(ks,1) ; 
    bw = bwconncomp(squeeze(ks(i,:,:))) ; plist = bw.PixelIdxList ; 
    maxsize = 1 ; maxind = 1 ; 
    for p=1:length(plist) ; if length(plist{p}) > maxsize ; maxind = p ; maxsize = length(plist{p}) ; end ; end 
    slice = zeros(size(ks,2),size(ks,3)) ; 
    if ~isempty(plist) ; slice(bw.PixelIdxList{maxind}) = 1 ; end
    remks(i,:,:) = slice ; 
end 
 %remks = medfilt3(remks) ; remks = medfilt3(remks) ; remks = medfilt3(remks) ;  
 %}
 
 
%create the outer shells (6-7) put in dils
dilmin = imdilate(remks,strel(ones(17,17,17)))-imdilate(remks,strel(ones(11,11,11))) ; 
bw = bwconncomp(dilmin) ; outer = zeros(size(dilmin)) ; outer(bw.PixelIdxList{1}) = 1 ; 
dil_outer = imdilate(outer,strel(ones(5,5,5))) ; outer = imerode(outer,strel(ones(3,3,3))) ; 

notouter = outer==0 ; inside = imdilate(remks,strel(ones(3,3,3))) ; 
for i=1:20 ; inside = inside.*notouter ; inside = imdilate(inside,strel(ones(3,3,3))) ; end
current_outer = outer ; 
for i=1:10 ; 
    dilout = imdilate(current_outer,strel(ones(3,3,3))) ; 
    dilout = dilout.*inside ; 
    if i==1 ; dils(i,:,:,:) = dilout ; else dils(i,:,:,:) = (dilout-prevdil) > 0 ; end
    current_outer = dilout ; 
    prevdil = dilout ; 
end

%{
icount = 1 ; dilprev = zeros(size(simg)) ; clear dils ; 
for i=3:2:15
   dili = imdilate(remks,strel(ones(i,i,i))) ; 
   multdil = dili.*dil_outer ; 
   if icount == 1 
      dils(icount,:,:,:) = multdil ;  
   elseif icount > 1 ;
      dils(icount,:,:,:) = (multdil - dilprev) > 0 ;      
   end
   dilprev = multdil + dilprev ; 
   icount = icount + 1 ;
end
%}

gs = double(simg-imfilter(simg,fspecial('gaussian',21,21))) ; 

% create the electrode pancake at different distances from the scalp
clear vqs ; 
max_xp = 2.5 ; min_xp = -2.5 ; max_yp = 2.5 ; min_yp = -2.5 ; 
for i=1:10 ; indimg = squeeze(dils(i,:,:,:)) ; 
inds = find(indimg==1) ;  
[x,y,z] = ind2sub(size(indimg),inds) ; 
% find the center of mass of the skull image:
[cx,cy,cz] = centmass3(simg) ; 
xdiffs = x-cx ; ydiffs = y-cy ; zdiffs = z-cz ; 
[th,phi,r] = cart2sph(xdiffs,ydiffs,zdiffs) ; 
allth{i} = th ; allphi{i} = phi ; allr{i} = r ; 
thbrain = zeros(size(simg)) ; thbrain(inds) = th ; 
phibrain = zeros(size(simg)) ; phibrain(inds) = phi ; 
rbrain = zeros(size(simg)) ; rbrain(inds) = r ; 
[xp,yp] = pol2cart(th,max(phi)-phi) ; % theta and rho are the phi and theta brains respectively
nsteps = 600 ; xsteps = min_xp:(max_xp-min_xp)/nsteps:max_xp ; ysteps = min_yp:(max_yp-min_yp)/nsteps:max_yp ;
[xg,yg] = meshgrid(xsteps,ysteps) ; 
vq = griddata(xp,yp,double(gs(inds)),xg,yg) ; vq(isnan(vq)) = 0 ; 
vqs(i,:,:) = vq ; 
end

clear rgb ; expo = 2 ;
rgb(:,:,1) = squeeze(round(mat2gray(mean(vqs(5,:,:),1)).^expo*255)) ;
rgb(:,:,2) = squeeze(round(mat2gray(mean(vqs(4,:,:),1)).^expo*255)) ;
rgb(:,:,3) = squeeze(round(mat2gray(mean(vqs(3,:,:),1)).^expo*255)) ;
figure,imagesc(rgb/255) ;


% do stuff in java...
%%% get the EEG electrodes from eeglab
EEG = pop_loadset('C:\Vision\Raw Files\noepi.set') ; 
rlabs = {EEG.chanlocs.labels} ; 
% get the electrode coordinates and perform the inverse transform on them
% to yield MRI space electrode locations
[labs,coords] = load_elocs ; 
for i=1:size(labs,2) ; labinds(i) = find(strcmpi(labs{i},rlabs)) ; end
rcopy = zeros(size(rgb)) ; m2g = zeros(size(rcopy,1),size(rcopy,2)) ; 
for i=1:size(coords,1)
    rcopy(coords(i,2)-1:coords(i,2)+1,coords(i,1)-1:coords(i,1)+1,:) = labinds(i) ; 
    m2g(coords(i,2)-1:coords(i,2)+1,coords(i,1)-1:coords(i,1)+1) = labinds(i) ; 
end
elecbrain = zeros(size(simg)) ; 
for elec=1:length(labinds) 
scalpinds = find(m2g==labinds(elec)) ; 
% pol2cart(theta,r)
xs = xg(scalpinds) ; ys = yg(scalpinds) ; 
[inv_theta,inv_phi] = cart2pol(xs,ys) ; inv_phi = 1.5708-inv_phi ; 
clear lowinds
for i=1:size(inv_theta,1)
    theta_i = inv_theta(i) ;  phi_i = inv_phi(i) ; 
    clear minjs mininds % find the closest theta,phi coordinate to the current scalp point 
    for j=1:size(allth,2)
        thdiffs = theta_i - allth{j} ; 
        phidiffs = phi_i - allphi{j} ; 
        sqrdiffs = sqrt(thdiffs.^2 + phidiffs.^2) ; 
        mis = find(sqrdiffs==min(sqrdiffs)) ;  mininds(j) = mis(1) ; 
        minjs(j) = sqrdiffs(mininds(j)) ; 
    end
    lowest_j = find(minjs==min(minjs)) ; lowest_j = lowest_j(1) ;  
    lowest_j_ind = mininds(lowest_j) ; 
    lowinds(i,:) = [lowest_j,lowest_j_ind] ; 
end
clear sph_pts
for i=1:size(lowinds,1) ; sph_pts(i,:) = [allth{lowinds(i,1)}(lowinds(i,2)),allphi{lowinds(i,1)}(lowinds(i,2)),allr{lowinds(i,1)}(lowinds(i,2))] ; end
[xd,yd,zd] = sph2cart(sph_pts(:,1),sph_pts(:,2),sph_pts(:,3)) ; 
coordx = uint8(xd+cx) ; coordy = uint8(yd+cy) ; coordz = uint8(zd+cz) ; 
for i=1:size(coordx,1) ; elecbrain(coordx(i),coordy(i),coordz(i)) = labinds(elec) ; end 
disp(elec) ;
end

s.img = elecbrain ; 
save_untouch_nii(s,'elecbrain.nii.gz') ;

head = load_untouch_nii('fast_t1/fast_seg_1.nii.gz') ; 
for i=1:size(labinds,2) ; [ecx(labinds(i)),ecy(labinds(i)),ecz(labinds(i))] = centmass3(elecbrain==labinds(i)) ; end
cmass{1} = ecx ; cmass{2} = ecy ; cmass{3} = ecz ; 

smoothmod = smoothn(modbrain) ; 
modgm = double(head.img).*smoothmod ; 
subplot(1,2,1) ; imagesc(squeeze(modgm(:,:,45)).*squeeze(simg(:,:,45))) ; 
subplot(1,2,2) ; imagesc(rot90(squeeze(modgm(110,:,:)).*squeeze(simg(110,:,:)))) ; 

modbrain = zeros(size(elecbrain)) ;
for i=1:size(labinds,2) ; modbrain(cmass{1}(labinds(i)),cmass{2}(labinds(i)),cmass{3}(labinds(i))) = modulation(labinds(i)) ; end

vcom = [120,200,35] ; diffx = cmass{1}-vcom(1) ;diffy = cmass{2}-vcom(2) ; diffz = cmass{3}-vcom(3) ; 
sumdiffs = sqrt(diffx.^2+diffy.^2+diffz.^2) ; 

plot(sumdiffs,modulation,'o','LineWidth',3) ; xlabel('electrode distance(mm) from V1') ; ylabel('SSVEP modulation (db)');





