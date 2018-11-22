clear all ; close all ; 

cd c:/shared/glaucoma/laurent/func/DICOM ; ls
times = load('c:\shared\glaucoma\laurent\laurent_stim_right.mat') ; 
times = times.stimtimes ; 
delay = load_untouch_nii('delay_in_1.nii.gz') ; 
delayimg = delay.img(:,:,:,1) ; 
% load the registered correlation brain => threshold, dilate and
% clusterize.
corrbrain = load_untouch_nii('corrbrain_1.nii.gz') ; 
dilcorrs = (imdilate(corrbrain.img>.2,strel(ones(3,3,3)))) ; 
zcorrs = zeros(size(dilcorrs)) ; 
bw = bwconncomp(dilcorrs) ; 
pix = bw.PixelIdxList ; 
lengths = cellfun(@length,pix) ; largest = find(lengths==max(lengths)) ; 
zcorrs(pix{largest}) = 1 ; 

% stimulus conditions:
outer = load('c:/shared/glaucoma/outer') ; outer = outer.outer ; 

startt = times(1,:) ; 
fmri = load_untouch_nii('bp_2.nii.gz') ; 
fimg = fmri.img ; 
resfimg = zeros(size(fimg,1),size(fimg,2),size(fimg,3),size(fimg,4)*2) ; 
for i=1:size(fimg,1) ; disp(i) ; 
    for j=1:size(fimg,2)
        for k=1:size(fimg,3)
            tsijk = squeeze(fimg(i,j,k,:)) ; 
            resfimg(i,j,k,:) = imresize(tsijk,[820,1]) ; 
        end
    end
end

TR = 1.24/2 ; 
startTR = round(startt./TR) ; 

% stimulus triggers 1=right wedge, moving CW, 2=annulus moving inwards,
% 3=right wedge, moving CCW, 4= annulus, moving outwards.
hrf = spm_hrf(1.24/2) ; % offset is 5 in the canonical HRF at this TR
clear epochs
inds = repmat(1:4,[1,4]) ; 
indcounts = ones(1,4) ; 
for i=1:length(startTR)   
    epochs(inds(i),indcounts(inds(i)),:,:,:,:) = squeeze(resfimg(:,:,:,startTR(i)-round(5/TR):startTR(i)+round(30/TR)+10)) ; 
    indcounts(inds(i)) = indcounts(inds(i)) + 1 ;   
end
mepochs = squeeze(mean(epochs,2)) ; 

zinds = find(zcorrs>=0)  ;
[mx,my,mz] = ind2sub(size(zcorrs),zinds) ; 


% FOR THE ANGLES
ntrs = 25/TR ; anglestep = 360/ntrs ; 
angles1 = mod(337.5:-anglestep:-22.5,360) ; 
angles2 = mod(-22.5:anglestep:337.5,360) ; 
allangles = unique([angles1,angles2]) ; 
anglebrain = zeros(size(epochs,3),size(epochs,4),size(epochs,5),length(allangles)) ; 
ts = zeros(length(mx),length(allangles)) ; 
for i=1:length(angles1) ; disp(angles1(i)) ; 
    indexi = find(allangles==angles1(i)) ;
    for j=1:length(zinds)
        delayj = delayimg(mx(j),my(j),mz(j)) ;
        if delayj > 7 ; delayj = 7 ; end
        oset = round(delayj./TR)+round(5/TR) ; 
        voxj = squeeze(mepochs(1,mx(j),my(j),mz(j),:)) ; 
        ts(j,indexi) = ts(j,indexi) + (voxj(i+oset)) ;  % - mean(voxj(5:10))) ; 
    end
end
for i=1:length(angles2) ; disp(angles2(i)) ; 
     indexi = find(allangles==angles2(i)) ;
    for j=1:length(zinds)
        delayj = delayimg(mx(j),my(j),mz(j)) ;
        if delayj > 7 ; delayj = 7 ; end
        oset = round(delayj./TR)+round(5/TR) ; 
        voxj = squeeze(mepochs(3,mx(j),my(j),mz(j),:)) ; 
        ts(j,indexi) = ts(j,indexi) + (voxj(i+oset)) ; %- mean(voxj(5:10))) ; 
    end
end

fullangles = zeros(size(fimg,1),size(fimg,2),size(fimg,3),size(ts,2)) ; 
for i=1:size(ts,1)
    fullangles(mx(i),my(i),mz(i),:) = ts(i,:) ;    
end
newangles = fullangles(:,:,:,1:2:end-1) + fullangles(:,:,:,2:2:end) ; 
cat = load_untouch_nii('cat.nii.gz') ; 
cat.img = newangles ; 
save_untouch_nii(cat,'anglebrain_2.nii.gz') ; 
for i=1:size(newangles,1) ; for j=1:size(newangles,2) ; for k=1:size(newangles,3) ; maxrhos(i,j,k) = find(squeeze(newangles(i,j,k,:))==max(squeeze(newangles(i,j,k,:))),1) ; end ; end ; end
f1 = load_untouch_nii('f1_1.nii.gz') ; f1.img = maxrhos ; save_untouch_nii(f1,'maxangles_2.nii.gz') ; 



% FOR THE RADIUS
ntrs = 25/TR ; rhostep = round(length(outer)/ntrs) ; 
rhos1 = outer(length(outer):-rhostep:1) ; 
rhos2 = outer(1:rhostep:length(outer)) ; 
oset = round(5/TR) + round(3.5/TR) ; 
allrhos = unique([rhos1,rhos2]) ; 
rhobrain = zeros(size(epochs,3),size(epochs,4),size(epochs,5),length(allrhos)) ; 
for i=1:length(rhos1) ; 
    indexi = find(allrhos==rhos1(i)) ; 
    rhobrain(:,:,:,indexi) = rhobrain(:,:,:,indexi) + squeeze(mepochs(2,:,:,:,i+oset)) ;  %- squeeze(mepochs(2,:,:,:,3)) ;  
end
for i=1:length(rhos2) ; 
    indexi = find(allrhos==rhos2(i)) ; 
    rhobrain(:,:,:,indexi) = rhobrain(:,:,:,indexi) + squeeze(mepochs(4,:,:,:,i+oset)) ; %- squeeze(mepochs(4,:,:,:,3)) ; 
end

newrho = rhobrain(:,:,:,1:2:end) + rhobrain(:,:,:,2:2:end) ; 
cat.img = newrho(:,:,:,1:40) ; 
save_untouch_nii(cat,'rhobrain_2.nii.gz') ; 
for i=1:size(newrho,1) ; for j=1:size(newrho,2) ; for k=1:size(newrho,3) ; maxrhos(i,j,k) = find(squeeze(newrho(i,j,k,:))==max(squeeze(newrho(i,j,k,:))),1) ; end ; end ; end
f1 = load_untouch_nii('f1_1.nii.gz') ; f1.img = maxrhos ; save_untouch_nii(f1,'maxrhos_2.nii.gz') ; 

%%% make some tuning curves



zinds = find(corrbrain.img>.3)  ;
[mx,my,mz] = ind2sub(size(zcorrs),zinds) ; 
for i=1:length(zinds)
    rhovox(i,:) = newrho(mx(i),my(i),mz(i),:) ;
    anglevox(i,:) = newangles(mx(i),my(i),mz(i),:) ; 
end







