
cd c:/shared/raw ; subs=dir('MONG*') ; 
EEG = pop_loadbv('C:\shared\MONG_01_RB\','MONG_01_RB_FIX_BOX.vhdr') ;
EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 

for sb=1:length(subs)
% get the electrode coordinates and perform the inverse transform on them
% to yield MRI space electrode locations
[labs,coords] = load_elocs(['C:\shared\elocs\',allelecs(sb).name,'.txt']) ; 
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


end
