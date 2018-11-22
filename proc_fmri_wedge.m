cd c:/shared/retino_test/ ; ls 

nii = load_untouch_nii('Test_Russell_2015_10_22-2_WIP_EEG-fMRI_S1.5_MB3_SH4__SENSE_12_1.nii') ; nimg = nii.img ; 
cd stimtrigs ; ls 
trigs = load('scanChecker_2.mat') ; 
logfile = cell2mat(trigs.logfile) ; 
trigtimes = logfile(1:2:end) ; trigtypes = logfile(2:2:end) ; 

TR = nii.hdr.dime.pixdim(5) ; 
trigtimes_TR = round(trigtimes./TR) ; 
hrf = spm_hrf(TR) ;
maxhrf = find(hrf==max(hrf)) ; 
maxhrfinds = maxhrf-1:maxhrf+1 ; 

ntimes = 10/TR ; 
anglestep = 360/ntimes ; 
startangles = [180+45, 90+45, 0+45, 270+45] ;
clear allangles ; 
for i=1:length(startangles)
    allangles(i,:) = startangles(i):anglestep:startangles(i)+360 ; 
end
allangles = mod(allangles,360) ; 

% find the angles for each trigger type...
resangles = reshape(allangles,[1,size(allangles,1)*size(allangles,2)]) ; 
for i=1:size(allangles,1) ; for j=1:size(allangles,2) ; arrinds(i,j) = sub2ind(size(allangles),i,j) ; end ; end
resarrinds = reshape(arrinds,[1,size(allangles,1)*size(allangles,2)]) ; 
[sv,si] = sort(resangles,'descend') ; 
sortinds = resarrinds(si) ; 
[sx,sy] = ind2sub(size(allangles),sortinds) ; 

angleindices = trigtypes ;
hrflen = round(length(hrf)/2) ; clear anglevals
angleindexcounts = ones(1,4) ; 
for i=1:length(trigtypes)
    trigtr = trigtimes_TR(i) ;
    trigtype = trigtypes(i) ;
    acti = nimg(:,:,:,trigtr:trigtr+hrflen+ceil(ntimes)) ;
    basei = nimg(:,:,:,trigtr-1) ;
    angleindex = angleindices(i) ; 
    for j=1:size(allangles,2)
        anglevals(angleindex,angleindexcounts(angleindex),j,:,:,:) = squeeze(mean(acti(:,:,:,maxhrfinds+j),4))-squeeze(mean(basei,4)) ; 
    end
    angleindexcounts(angleindex) = angleindexcounts(angleindex) + 1 ; 
end
mangles = squeeze(mean(anglevals,2))./squeeze(std(anglevals,0,2)) ; mangles(isnan(mangles)) = 0 ; mangles(isinf(mangles)) = 0 ; 
resmangles = reshape(mangles,[1,numel(mangles)]) ; resmangles(abs(zscore(resmangles))>5) = 0 ; 
mangles = reshape(resmangles,size(mangles)) ; 
sortmangles = zeros(size(mangles,3),size(mangles,4),size(mangles,5),size(mangles,1)*size(mangles,2)) ; 
for i=1:length(sx)
    sortmangles(:,:,:,i) = squeeze(mangles(sx(i),sy(i),:,:,:)) ;
end
save_nii(make_nii(sortmangles),'sortmangles.nii.gz') ; 


