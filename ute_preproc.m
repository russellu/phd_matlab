%{
homedir = 'c:/shared/badger/felix/ute' ; 
cd(homedir) ; ls 
raw = load_untouch_nii('res_ute.nii.gz') ; 
gs = raw.img - imfilter(raw.img,fspecial('gaussian',41,21)) ; 
raw.img = gs ; save_untouch_nii(raw,'gs.nii.gz') ; 
%}

clear all ; close all ; disp('loading data...') ; 
%cd c:/shared/mong_05_SG/ ; ls ; s = load_untouch_nii('res_ute.nii.gz') ; 
%cd c:/shared/UTE ; ls ; s = load_untouch_nii('res_lin.nii.gz') ; 
cd c:/shared/badger/felix/ute ;
EEG = pop_loadbv('C:\shared\raw_eeg\MONG_01_RB','MONG_01_RB_FIX_BOX.vhdr') ;
EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
rlabs = {EEG.chanlocs.labels} ; 

for el=1%:length(allelecs) ; 
cd(['c:/shared/badger/felix/ute']) ; %s = load_untouch_nii('res_ute.nii.gz') ; 
elayers = load_untouch_nii('intensitylayers.nii.gz') ;  eimg = elayers.img ;   
for i=1:max(max(max(eimg))) ; dils(i,:,:,:) = eimg==i ; end
%t1 = load_untouch_nii('t1.nii.gz') ; t1img = double(t1.img - imfilter(t1.img,fspecial('gaussian',41,41))) ; 
s = load_untouch_nii('res_ute.nii.gz') ; 
simg = s.img ; gs = double((simg)-imfilter((simg),fspecial('gaussian',90,90))) ; 
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

for i=1:max(max(max(eimg))) ; 
    dqs(i,:,:) = (squeeze(vqs(i,:,:)))-imfilter((squeeze(vqs(i,:,:))),fspecial('gaussian',25,25)) ;  
end
expo = 1 ; 
dqs(:,155:170,220:240) = 0 ; 
rgbs(:,:,3) = uint8(mat2gray(squeeze(mean(dqs(5:6,:,:),1)).^expo)*255) ;  
rgbs(:,:,2) = uint8(mat2gray(squeeze(mean(dqs(3:4,:,:),1)).^expo)*255) ;  
rgbs(:,:,1) = uint8(mat2gray(squeeze(mean(dqs(1:2,:,:),1)).^expo)*255) ;  
figure,imagesc(rgbs)
allrgbs(el,:,:,:) = rgbs ; 
allvqs(el,:,:,:) = vqs ; 
imwrite(rgbs,['c:/shared/UTE_pngs/','felix_native','.png']) ;

% do stuff in java...

[labs,coords] = load_elocs(['C:\shared\elocs\felix_native.txt']) ; 
for i=1:size(labs,2)-2 ; labinds(i) = find(strcmpi(labs{i},rlabs)) ; end
rcopy = zeros(size(rgbs)) ; m2g = zeros(size(rcopy,1),size(rcopy,2)) ; 
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
        % rdiffs => need an r value here also...get the indices from theta
        % and phi and then grab the r-values from those indices?
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
coordx = round(xd+cx) ; coordy = round(yd+cy) ; coordz = round(zd+cz) ; 
for i=1:size(coordx,1) ; elecbrain(coordx(i),coordy(i),coordz(i)) = labinds(elec) ; end 
disp(elec) ;
end
%elecbrain = imdilate(elecbrain,strel(ones(3,3,3))) ; 
%s.img = elecbrain ;
%save_untouch_nii(s,'elecbrain_2.nii.gz') ;
%head = load_untouch_nii('fast_t1/fast_seg_1.nii.gz') ;  
cd(['c:/shared/badger/felix/ute'])
elayers.img = uint8(elecbrain>0) ; save_untouch_nii(elayers,'elecbrain.nii.gz')  ;
clear cmass ;
for i=1:size(labinds,2) ; [ecx(labinds(i)),ecy(labinds(i)),ecz(labinds(i))] = centmass3(elecbrain==labinds(i)) ; end
cmass{1} = ecx ; cmass{2} = ecy ; cmass{3} = ecz ; 
save('cmass','cmass') ; 
elabels = rlabs(labinds) ; save('elabels','elabels') ; 

% make the skin mask
nlayers = load_untouch_nii('native_layers.nii.gz') ; 
nlayersimg = nlayers.img ; 
skin = imdilate(nlayersimg==1,strel(ones(5,5,5))).*(nlayersimg==0) ;
nlayers.img = skin ; save_untouch_nii(nlayers,'skin.nii.gz') ; 
% use the skin to project vectors to the center of mass and pinpoint the
% electrode intersection with the surface of the head
[cx,cy,cz] = centmass3(skin) ; 
skin(cx-5:cx+5,cy-5:cy+5,cz-5:cz+5) = 1 ; 
emass = [ecx;ecy;ecz] ; 
diffx = ecx-cx ; diffy = ecy-cy ; diffz = ecz-cz ; 
diffvecs = [diffx;diffy;diffz] ;
unitdiffs = diffvecs./repmat(sqrt(sum(diffvecs.^2,1)),[3,1]) ; 
% make the lines radiating outwards from the center of mass to each
lsize = 150 ;
lines = repmat((1:lsize)',[1,3,64]) ; 
repunits = zeros(1,3,64) ; repunits(1,:,:) = unitdiffs ; repunits = repmat(repunits,[lsize,1,1]) ; 
multlines = lines.*repunits ; 
multlines(:,1,:) = multlines(:,1,:) + cx ; multlines(:,2,:) = multlines(:,2,:) + cy ; multlines(:,3,:) = multlines(:,3,:) + cz ; 
multlines = round(multlines) ; 
skinlines = zeros(size(skin)) ; cskinlines = zeros(size(skin)) ; 

elecinds = [1:31,33:64] ; 
for i=1:size(multlines,1) ; 
    for j=1:length(elecinds) ; 
        ix = multlines(i,1,elecinds(j)) ; iy = multlines(i,2,elecinds(j)) ; iz = multlines(i,3,elecinds(j)) ;
        if ix>0 & ix<size(skin,1) & iy>0 & iy<size(skin,2) & iz>0 & iz<size(skin,3)
            skinlines(ix,iy,iz) = 1 ; 
            cskinlines(ix,iy,iz) = (find(labinds==(elecinds(j)))) ; % assign the line a color based on the electrode you are projecting to
        end
    end ; 
end 
nlayers.img = skinlines ; save_untouch_nii(nlayers,'skinlines.nii.gz') ; 
skin(cx-5:cx+5,cy-5:cy+5,cz-5:cz+5) = 0 ; 

skinsurf = load_untouch_nii('skinsurf.nii.gz') ; 
skinter = skinlines.*skinsurf.img  ; cskinter = imdilate(cskinlines.*skin,strel(ones(7,7,7))) ; 

nlayers.img = imdilate(skinter,strel(ones(3,3,3))) ; save_untouch_nii(nlayers,'skinter.nii.gz') ; 
save('elabels','elabels') ; 
nlayers.img = cskinter ; save_untouch_nii(nlayers,'cskinter.nii.gz') ; 
%}
end
