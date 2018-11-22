clear all ; close all ; 
sujets = {'sujet1','sujet2','sujet3','sujet4','sujet5','sujet6','sujet8','sujet9','sujet10','sujet11','sujet12'} ; 
coords = {[-54.247, -34.634, 100.716],[-56.828 -73.772 29.952],[ -57.938 -46.949 61.904],[-54.003 -40.933 57.548],[-62.218 -71.093 14.529],[-51.816	-72.589	72.933],[-56.472	-56.854	53.331],[-53.761 -49.228 46.651]...
    [-50.379 -70.730 56.868],[-56.330 -60.345 53.377],[-58.909 -46.796 32.446]} ; 

for s=1:length(sujets) ; 
cd(['c:/shared/claudie/',sujets{s}]) ; ls 
mni = load_untouch_nii('deob_t1.nii.gz') ; 

xyz = coords{s} ; 
offx = mni.hdr.hist.qoffset_x ; % right to middle
offy = mni.hdr.hist.qoffset_y ; % middle voxel to back of image
offz = mni.hdr.hist.qoffset_z ; 
disp(['offx = ',num2str(offx), ' offy = ',num2str(offy), ' offz = ',num2str(offz)]) ; 
voxcoords = xyz ;
voxcoords(:,1) = -(voxcoords(:,1)-offx) ; 
voxcoords(:,2) = -(voxcoords(:,2)-offy) ; 
voxcoords(:,3) = voxcoords(:,3) -offz ; 
voxcoords = round(voxcoords') ; 
elecs = zeros(size(mni.img)) ; 

for i=1:size(voxcoords,2)
   elecs(voxcoords(1,i),voxcoords(2,i),voxcoords(3,i)) = i ;    
end

elecs = imdilate(elecs,strel(ones(3,3,3))) ; 
mni.img = elecs ; 
save_untouch_nii(mni,'parietal_zone.nii.gz') ; 
end

