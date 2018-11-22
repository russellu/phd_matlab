clear all ; close all ; disp('loading data...') ; 
%cd c:/shared/mong_05_SG/ ; ls ; s = load_untouch_nii('res_ute.nii.gz') ; 
%cd c:/shared/UTE ; ls ; s = load_untouch_nii('res_lin.nii.gz') ; 
cd c:/shared/raw ; allelecs = dir('*MONG*') ; 

for el=2:length(allelecs) ; 
    
cd(['c:/shared/raw/',allelecs(el).name]) ; s = load_untouch_nii('elecs.nii.gz') ; 
%t1 = load_untouch_nii('t1_in_elecs.nii.gz') ; t1im = double(t1.img)-imfilter(double(t1.img),fspecial('gaussian',90,90))  ; 
padamt = 100 ; 
simg = double(pad3d(double(s.img),padamt)) ; 
gs = (sqrt(simg)-imfilter(sqrt(simg),fspecial('gaussian',90,90))) ; 
[cx,cy,cz] = centmass3(simg) ;
% find the enclosing sphere (smallest)
disp('creating grid') ; 
[xn,yn,zn] = ndgrid(-cx:size(simg,1)-cx-1,-cy:size(simg,2)-cy-1,-cz:size(simg,3)-cz-1) ; 
disp('creating bounding sphere...') ; 
minrad = 120 ; 
minsphere = sqrt(xn.^2+yn.^2+zn.^2) > minrad-1 & sqrt(xn.^2+yn.^2+zn.^2) < minrad+1 ; 

disp('finding line indices...') ; 
% robust scalp detection 
cmass = [cx,cy,cz] ; 
clear brainlines ; 
[xi,yi,zi] = ind2sub(size(minsphere),find(minsphere==1)) ;    
% the difference in voxels from each boundary point to the centmass
massdiffs = [xi,yi,zi] - repmat([cx,cy,cz],[size([xi,yi,zi],1),1]) ;
% the unit vector between each boundary point and the centmass
normdiffs = massdiffs./sqrt(repmat(sum(massdiffs.^2,2),[1,3])) ; 
% use the unit vector in increments to go from center of mass to
linelength = 150 ; 
cs_x = ceil(repmat(cumsum(ones(1,linelength)),[size(normdiffs,1),1]).*repmat(normdiffs(:,1),[1,linelength]) + cx) ; 
cs_y = ceil(repmat(cumsum(ones(1,linelength)),[size(normdiffs,1),1]).*repmat(normdiffs(:,2),[1,linelength]) + cy) ; 
cs_z = ceil(repmat(cumsum(ones(1,linelength)),[size(normdiffs,1),1]).*repmat(normdiffs(:,3),[1,linelength]) + cz) ; 
lineinds = sub2ind(size(minsphere),cs_x,cs_y,cs_z) ; 
resinds = reshape(lineinds,[1,size(lineinds,1)*size(lineinds,2)]) ; 
resvals = gs(resinds) ; 
brainlines = reshape(resvals,size(lineinds)) ; 
maxes = max(brainlines,[],2) ; 
disp('finding maximum indices...') ; 
for i=1:size(brainlines,1) ; maxinds(i) = lineinds(i,find(brainlines(i,:)==maxes(i),1)) ; end
maxbrain = zeros(size(gs)) ;
for i=1:length(maxinds) ; maxbrain(maxinds(i)) = maxbrain(maxinds(i))+1 ; end
newskull = maxbrain(padamt:size(maxbrain,1)-padamt-1,padamt:size(maxbrain,2)-padamt-1,padamt:size(maxbrain,3)-padamt-1) ; 
newskull = imdilate(newskull>0,strel(ones(7,7,7))) ; ccomps = bwconncomp(newskull) ; skullinds = ccomps.PixelIdxList{1} ; 
zskull = zeros(size(newskull)) ; zskull(skullinds) = 1 ; 
zskull = (imfilter(zskull,fspecial('gaussian',7,3))>0) ; 
bw2 = imfill(zskull,[1,1,1],6) ; 
outer = bw2+zskull*2 == 1 ;  erodeout = imerode(outer,strel(ones(9,9,9))) ;
layers = zeros(size(erodeout)) ; 
for i=1:10 ; 
    newerode = imdilate(erodeout,strel(ones(3,3,3))) ; 
    erodediff = (erodeout - newerode) == -1 ;  
    layers = layers + erodediff*i ; 
    dils(i,:,:,:) = erodediff ; 
    erodeout = erodeout + erodediff ; 
end
s.img = layers ; save_untouch_nii(s,'layers.nii.gz') ; 


s = load_untouch_nii('elecs.nii.gz') ; 
simg = s.img ; gs = double((simg)-imfilter((simg),fspecial('gaussian',50,50))) ; 
s.img = gs ; save_untouch_nii(s,'gs.nii.gz') ; 
clear vqs ; 
max_xp = 2.5 ; min_xp = -2.5 ; max_yp = 2.5 ; min_yp = -2.5 ; 
for i=1:size(dils,1) ; indimg = squeeze(dils(i,:,:,:)) ; 
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
nsteps = 550 ; xsteps = min_xp:(max_xp-min_xp)/nsteps:max_xp ; ysteps = min_yp:(max_yp-min_yp)/nsteps:max_yp ;
[xg,yg] = meshgrid(xsteps,ysteps) ; 
vq = griddata(double(xp),double(yp),gs(inds),double(xg),double(yg)) ; vq(isnan(vq)) = 0 ; 
vqs(i,:,:) = vq ; 
end
%figure,
%for i=1:10 ; subplot(3,4,i) ;imagesc(squeeze(vqs(i,:,:))) ; colormap gray ; end
expo = 1 ; 
rgbs(:,:,1) = uint8(mat2gray(squeeze(mean(vqs(6:7,:,:),1)).^expo)*255) ;  
rgbs(:,:,2) = uint8(mat2gray(squeeze(mean(vqs(4:5,:,:),1)).^expo)*255) ;  
rgbs(:,:,3) = uint8(mat2gray(squeeze(mean(vqs(2:3,:,:),1)).^expo)*255) ;  
allrgbs(el,:,:,:) = rgbs ; 
allvqs(el,:,:,:) = vqs ; 
 


   % figure, imagesc(drgbs) ; 
   % imwrite(drgbs,['c:/shared/UTE_pngs/',allelecs(el).name,'.png']) ; 
end


% do stuff in java...
%%% get the EEG electrodes from eeglab
EEG = pop_loadbv('C:\shared\MONG_01_RB\','MONG_01_RB_FIX_BOX.vhdr') ;
EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
rlabs = {EEG.chanlocs.labels} ; 
% get the electrode coordinates and perform the inverse transform on them
% to yield MRI space electrode locations
[labs,coords] = load_elocs(['C:\shared\elocs\',allelecs(el).name,'.txt']) ; 
for i=1:size(labs,2)-2 ; labinds(i) = find(strcmpi(labs{i},rlabs)) ; end
rcopy = zeros(size(drgbs)) ; m2g = zeros(size(rcopy,1),size(rcopy,2)) ; 
for i=1:size(coords,1)-2
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
s.img = elecbrain ; cd(['c:/shared/raw/',allelecs(el).name])
save_untouch_nii(s,'elecbrain.nii.gz') ;
%head = load_untouch_nii('fast_t1/fast_seg_1.nii.gz') ; 
clear cmass ;
for i=1:size(labinds,2) ; [ecx(labinds(i)),ecy(labinds(i)),ecz(labinds(i))] = centmass3(elecbrain==labinds(i)) ; end
cmass{1} = ecx ; cmass{2} = ecy ; cmass{3} = ecz ; 
save('cmass','cmass') ; 

%{
smoothmod = smoothn(modbrain) ; 
modgm = double(head.img).*smoothmod ; 
subplot(1,2,1) ; imagesc(squeeze(modgm(:,:,45)).*squeeze(simg(:,:,45))) ; 
subplot(1,2,2) ; imagesc(rot90(squeeze(modgm(110,:,:)).*squeeze(simg(110,:,:)))) ; 
modbrain = zeros(size(elecbrain)) ;
for i=1:size(labinds,2) ; modbrain(cmass{1}(labinds(i)),cmass{2}(labinds(i)),cmass{3}(labinds(i))) = modulation(labinds(i)) ; end
vcom = [120,200,35] ; diffx = cmass{1}-vcom(1) ;diffy = cmass{2}-vcom(2) ; diffz = cmass{3}-vcom(3) ; 
sumdiffs = sqrt(diffx.^2+diffy.^2+diffz.^2) ; 
plot(sumdiffs,modulation,'o','LineWidth',3) ; xlabel('electrode distance(mm) from V1') ; ylabel('SSVEP modulation (db)');
%}