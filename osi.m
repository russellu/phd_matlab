cd c:/shared/badger_mri/russell/nii ; ls 
allangles = load_untouch_nii('allangles_t1.nii.gz') ; 
angimg = allangles.img ; 
maxangles = max(angimg,[],4) ; minangles = min(angimg,[],4) ; 
x = maxangles-minangles ; 
maxinds = zeros(size(maxangles)) ; 
for i=1:size(angimg,1)
    for j=1:size(angimg,2)
        for k=1:size(angimg,3)
            if maxangles(i,j,k) ~= 0
                maxinds(i,j,k) = find(squeeze(angimg(i,j,k,:))==squeeze(maxangles(i,j,k)),1) ; 
            end
        end
    end
end

v1clust = load_untouch_nii('v1_cluster.nii.gz') ; 

segimg = zeros(size(maxinds,1),size(maxinds,2),size(maxinds,3),20) ; 
icount = 1 ; 
for i=0:4.5:90
    segimg(:,:,:,icount) = (maxinds>i & maxinds<i+4.5).*v1clust.img ; 
    icount = icount + 1 ; 
end
save_nii(make_nii(segimg),'segimg.nii.gz') ; 



% compare to the vasculature 
ved = load_untouch_nii('veins.nii.gz') ; 
nonz = find(x>30) ; 
nonv = find(ved.img>5000) ; 
goodinds = intersect(nonz,nonv) ; 
vedvals = ved.img(goodinds) ; 
xvals = x(goodinds) ; 





















