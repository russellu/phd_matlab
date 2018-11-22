clear all ; close all ; 
cd c:/brainstorm3/
ls 
xyz = load('xyz_mni.mat') ; xyz = xyz.xyz_mni; 
labs = load('labs.mat') ; labs = labs.labs ; 
cd c:/shared/ATLASES ; ls 
mni = load_untouch_nii('MNI152_T1_1mm.nii.gz') ; 

offx = mni.hdr.hist.qoffset_x ; % right to middle
offy = mni.hdr.hist.qoffset_y ; % middle voxel to back of image
offz = mni.hdr.hist.qoffset_z ; 

xyz = xyz*1000 ; 

cx = xyz(1,1) ; 
cy = xyz(2,1) ; 
cz = xyz(3,1) ; 

voxcoords = xyz ;
voxcoords(1,:) = -(voxcoords(1,:) - offx) ; 
voxcoords(2,:) = (voxcoords(2,:) - offy) ; 
voxcoords(3,:) = voxcoords(3,:) -offz ; 

rvc = round(voxcoords) ; 
elecs = zeros(size(mni.img)) ; 
for i=1:size(rvc,2)
   elecs(rvc(1,i),rvc(2,i),rvc(3,i)) = i ;    
end

elecs = imdilate(elecs,strel(ones(3,3,3))) ; 
mni.img = elecs ; 
save_untouch_nii(mni,'mni_elecs.nii.gz') ; 



cd c:/shared/ATLASES ; ls 
mask = load_untouch_nii('mnimask.nii.gz') ; 












