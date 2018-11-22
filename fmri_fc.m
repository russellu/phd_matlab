cd c:/shared/allfmris/sub_russell ; 
r = load_untouch_nii('common_1.nii.gz') ; 
rim = squeeze(mean(r.img,4)) ; 

brim = rim>100 ; 
binds = find(brim==1) ; 
[ix,iy,iz] = ind2sub(size(brim),binds) ; 

rimg = r.img ; 
pixt = zeros(length(ix),size(rimg,4)) ; 
for i=1:length(ix) ; pixt(i,:) = squeeze(rimg(ix(i),iy(i),iz(i),:)) ; end







