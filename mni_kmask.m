cd c:/shared/ATLASES
ls
a = load_untouch_nii('medgm_3mm.nii.gz') ; 

inds = find(a.img==1) ; 
[gx,gy,gz] = ndgrid(-size(a.img,1)/2:size(a.img,1)/2-1,-size(a.img,2)/2:size(a.img,2)/2-1,-size(a.img,3)/2:size(a.img,3)/2-1) ; 
for i=1:length(inds)
   indcoords(i,:) = [gx(inds(i)),gy(inds(i)),gz(inds(i))] ;  
end

kvals = kmeans(indcoords,200) ; 
kbrain = zeros(size(gx)) ; 
for i=1:length(kvals)
    kbrain(inds(i)) = kvals(i) ;   
end

a.img = kbrain ; save_untouch_nii(a,'kbrain.nii.gz') ; 







