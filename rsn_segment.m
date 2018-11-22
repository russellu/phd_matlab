cd C:\shared\epireg\rsns
nimgs = 0; 
all=dir('*gz'); 
for i=1:length(all)
    mri = load_untouch_nii(all(i).name); 
    epis{i} = mri.img; 
    nimgs = nimgs + size(mri.img,4); 
end

allimgs = zeros(size(mri.img,1),size(mri.img,2),size(mri.img,3),nimgs); 
imgcount = 1;
for i=1:length(epis)
    for j=1:size(epis{i},4)
        allimgs(:,:,:,imgcount) = epis{i}(:,:,:,j); 
        imgcount = imgcount + 1 ; 
    end
end

cd C:\shared\epireg\
meanimg = load_untouch_nii('mean.nii.gz'); 
voxels = find(meanimg.img>600); 
[vx,vy,vz] = ind2sub(size(meanimg.img),voxels); 
ts = zeros(size(allimgs,4),length(vx)); 
for i=1:size(allimgs,4)
    for j=1:length(vx)
        ts(i,j) = allimgs(vx(j),vy(j),vz(j),i); 
    end
end

corrs = corr(ts') ;
[sv,si] = sort(corrs,2,'descend'); 
corrimgs = zeros(size(allimgs)); 
for i=1:size(si,1)
    meani = squeeze(mean(allimgs(:,:,:,si(i,2:12)),4)); 
    corrimgs(:,:,:,i) = meani; 
end

allimgs = corrimgs; 
ts = zeros(size(allimgs,4),length(vx)); 
for i=1:size(allimgs,4)
    for j=1:length(vx)
        ts(i,j) = allimgs(vx(j),vy(j),vz(j),i); 
    end
end
corrs = corr(ts') ;
[sv,si] = sort(corrs,2,'descend'); 
corrimgs = zeros(size(allimgs)); 
for i=1:size(si,1)
    meani = squeeze(mean(allimgs(:,:,:,si(i,2:12)),4)); 
    corrimgs(:,:,:,i) = meani; 
end

allimgs = corrimgs; 
ts = zeros(size(allimgs,4),length(vx)); 
for i=1:size(allimgs,4)
    for j=1:length(vx)
        ts(i,j) = allimgs(vx(j),vy(j),vz(j),i); 
    end
end
corrs = corr(ts') ;
[sv,si] = sort(corrs,2,'descend'); 
corrimgs = zeros(size(allimgs)); 
for i=1:size(si,1)
    meani = squeeze(mean(allimgs(:,:,:,si(i,2:12)),4)); 
    corrimgs(:,:,:,i) = meani; 
end


cd C:\shared\epireg\rsns
template = load_untouch_nii('allsubs.nii.gz'); 
template.img = corrimgs; 
save_untouch_nii(template,'corrimgs.nii.gz'); 

