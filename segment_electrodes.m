%%% segment EEG electrodes based off of UTE image intensity and the
%%% hand-labeled points from "hand_labeling.m"
%%% outputs final electrode coordinates based off of the segmentation
%%% procedure, in UTE coordinates
% author Russell Butler russell.buttler@usherbrooke.ca
clear all ; close all ; 

handLabel_output_path = 'c:/shared/ute_output/russ' ; % should be the same as path_output from hand_labeling.m
cd(handLabel_output_path) ; 
disp('loading raw data...') ; 
rute = load_untouch_nii('res_ute.nii.gz') ; ruteorig = double(rute.img) ;  
fmask = load_untouch_nii('finalmask.nii.gz') ; maskimg = imdilate(fmask.img>0,strel(ones(3,3,3))) ; maskimg = ~maskimg ; 
outerint = maskimg.*ruteorig ; 

boxl = 10 ; % size of half the box length (20x20x20mm)
% pad the output images, so when we take the electrode's neighbourhood it
% doesn't go out of bounds
outerint = pad3d(outerint,boxl) ; 
boximg = zeros(size(outerint)) ; 
squareimg = zeros(size(outerint)) ; % to visualize the boxes
pointimg = zeros(size(outerint)) ; % to visualize the hand-labeled points
segimg = zeros(size(outerint)) ; 
segcoordimg = zeros(size(outerint)) ; 
rndcoordimg = zeros(size(outerint)) ; 
bincentavg = zeros(size(outerint)) ; 
colorcentavg = zeros(size(outerint)) ; 



for c=1:30 % for each randomized coordinate set of 65 electrodes
    
coords = load(['mricoords_',num2str(c),'.mat']) ; % get the coordinates
if c==1
   origcoords = coords.mricoords + boxl ; 
end

coords = coords.mricoords + boxl ; % add the pad offset to each coordinate

for i=1:size(coords,2) % for all electrodes
    boxi = outerint(coords(1,i)-boxl:coords(1,i)+boxl,coords(2,i)-boxl:coords(2,i)+boxl,coords(3,i)-boxl:coords(3,i)+boxl) ; % get the box
    if c==1 % if the original hand-labeled points, save the info as a volume for visualization later
        squareimg(coords(1,i)-boxl:coords(1,i)+boxl,coords(2,i)-boxl:coords(2,i)+boxl,coords(3,i)-boxl:coords(3,i)+boxl) = 1 ; 
        pointimg(coords(1,i),coords(2,i),coords(3,i)) = 1 ; 
    end
    [cx,cy,cz] = centmass3(boxi) ; 
    normcent(c,i,:) = [cx+coords(1,i),cy+coords(2,i),cz+coords(3,i)] ; % normcent is the center of mass without segmentation (ie, based solely on image intensity within the cube)
    resboxi = reshape(boxi,[1,numel(boxi)]) ; 
    [sv,si] = sort(resboxi,'descend') ; 
    resboxi(si(1:200)) = 1  ; resboxi(si(101:end)) = 0 ; % set the voxels to 0 or 1
    boxi = reshape(resboxi,size(boxi)) ; 
    [cx2,cy2,cz2] = centmass3(boxi) ; 
    bincent(c,i,:) = [cx2+coords(1,i),cy2+coords(2,i),cz2+coords(3,i)] ; % bincent is the center of mass WITH segmentation (ie, based on the binary image intensity after taking top 100 voxels)
    allboxes(i,:,:,:) = boxi ; 
    if c==1 % save the segmented coordinates and electrode clusters for visualization
        segimg(coords(1,i)-boxl:coords(1,i)+boxl,coords(2,i)-boxl:coords(2,i)+boxl,coords(3,i)-boxl:coords(3,i)+boxl) = boxi ; 
        segcoordimg(cx2+coords(1,i)-boxl,cy2+coords(2,i)-boxl,cz2+coords(3,i)-boxl) = 1 ; 
    end
    rndcoordimg(coords(1,i),coords(2,i),coords(3,i)) = 1 ; 
end
end

% create an image with all the centers of mass
boximg2 = zeros(size(boximg)) ; 
for i=1:size(normcent,1)
    for j=1:size(normcent,2)
        boximg(normcent(i,j,1)-boxl,normcent(i,j,2)-boxl,normcent(i,j,3)-boxl) = j ; 
        boximg2(bincent(i,j,1)-boxl,bincent(i,j,2)-boxl,bincent(i,j,3)-boxl) = j ; 
    end
end

meanbincent = round(squeeze(mean(bincent,1)))-boxl ; 
meanstds = mean(squeeze(std(bincent,0,1)),2) ; % mean std in a single subject across 30 permutations
for i=1:size(meanbincent,1)
    if meanstds(i) < 1
        bincentavg(meanbincent(i,1),meanbincent(i,2),meanbincent(i,3)) = 1 ; 
        colorcentavg(meanbincent(i,1),meanbincent(i,2),meanbincent(i,3)) = i ; 
    else
        bincentavg(origcoords(1,i),origcoords(2,i),origcoords(3,i)) = 1 ; 
        colorcentavg(origcoords(1,i),origcoords(2,i),origcoords(3,i)) = i ; 
    end
end

% undo the zero pad, so it matches the dimensions of the original nifti
boximg = boximg(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 
boximg2 = boximg2(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 
squareimg = squareimg(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 
pointimg = pointimg(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 
segimg = segimg(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 
segcoordimg = segcoordimg(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 
rndcoordimg = rndcoordimg(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 
bincentavg = bincentavg(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 
colorcentavg = colorcentavg(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 

% a bunch of output images, mostly for visualization purposes for the
% article figures
rute.img = boximg2>0 ; save_untouch_nii(rute,'boximg2.nii.gz') ; 
rute.img = boximg>0 ; save_untouch_nii(rute,'boximg.nii.gz') ; 
rute.img = squareimg ; save_untouch_nii(rute,'squareimg.nii.gz') ; 
rute.img = imdilate(pointimg,strel(ones(3,3,3))) ; save_untouch_nii(rute,'pointimg.nii.gz') ; 
rute.img = segimg ; save_untouch_nii(rute,'segimg.nii.gz') ; 
rute.img = imdilate(segcoordimg,strel(ones(3,3,3))) ; save_untouch_nii(rute,'segcoordimg.nii.gz') ; 
rute.img = rndcoordimg ; save_untouch_nii(rute,'rndcoordimg.nii.gz') ; 
rute.img = imdilate(bincentavg,strel(ones(3,3,3))) ; save_untouch_nii(rute,'bincentavg.nii.gz') ; 
rute.img = colorcentavg ; save_untouch_nii(rute,'colorcentavg.nii.gz') ; 

