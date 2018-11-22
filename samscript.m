clear all ; close all ; 
cd C:\shared\for_russell ; ls 

% get the 3 datasets (contra, rest, and ipsi)
contra_stim = load_untouch_nii('tof_R_leftBranch.nii') ; 
rest = load_untouch_nii('tof_rest_leftBranch.nii') ; 
ipsi = load_untouch_nii('tof_L_leftBranch.nii') ; 

% dilate the 3 images, to improve the intersection
dilcontra = imdilate(contra_stim.img,strel(ones(3,3,3))) ; 
dilrest = imdilate(rest.img,strel(ones(3,3,3))) ; 
dilipsi = imdilate(ipsi.img,strel(ones(3,3,3))) ; 

% visualize the intersections with/without dilations using MIPs
% there is some overlap, but in general without using imdilate there is very little
% overlap between the 3 images. 
figure,
subplot(2,1,1) ; imagesc((max(contra_stim.img,[],3)>0) + (max(rest.img,[],3)>0) + (max(ipsi.img,[],3)>0)) ; colorbar ; 
title('intersection of MIP of all 3 images, RED = all 3 overlap') ; 
subplot(2,1,2) ; imagesc((max(dilcontra,[],3)>0) + (max(dilrest,[],3)>0) + (max(dilipsi,[],3)>0)) ; colorbar ; 
title('intersection of MIP of all 3 dilated images, RED = all 3 overlap') ; 
suptitle('FIGURE 1') ; 

% get the indices where all 3 DILATED images overlap
% to do this, you need to threshold each image separately, and find the
% indices where the sum of all 3 thresholded images is equal to 3
inds = find(((dilcontra>0) + (dilrest>0) + (dilipsi>0)) == 3) ; 

% now, you can use these indices to subtract diameters from voxels where
% you know the 3 images overlap
contra_diams = dilcontra(inds) ; % get the contralateral diameters
rest_diams = dilrest(inds) ; % get the rest diameters
ipsi_diams = dilipsi(inds) ; % get the ipsi diameters
contra_minus_rest = contra_diams - rest_diams ; % subtraction

% visualize a histogram of the subtracted diameters (50 bins) 
figure,
subplot(2,1,1) ; hist(contra_minus_rest,50) ; xlabel('contra-rest (mm)') ; 
title('contra minus rest diameter differences') ; 
zimg = zeros(size(dilcontra)) ; zimg(inds) = contra_minus_rest ; 
subplot(2,1,2) ; imagesc(sum(zimg,3)) ; colorbar ; 
title('diameter differences contra minus rest (arbitrary units)')
suptitle('FIGURE 2') ; 

% get the indices where there is a dilation from rest to contra. you can
% change '0' to whatever you want, if you want to set a dilation threshold
dilation_indices = find(contra_minus_rest > 0) ; 

% get the diameters in the voxels where contralateral was greater than
% rest, and subtract rest-ipsi in the same voxels
restdiams_where_dilated = rest_diams(dilation_indices) ; % rest diams
ipsidiams_where_dilated = ipsi_diams(dilation_indices) ; % ipsi diams
rest_minus_ipsi = restdiams_where_dilated - ipsidiams_where_dilated ; % differences

% display the final result
figure,
hist(rest_minus_ipsi) ; xlabel('ipsi diameter change, rest-stim') ; 
[h,p,ci,stats] = ttest(rest_minus_ipsi) ; % t-test
title(['histogram of diameter changes [p=',num2str(p),', t=',num2str(stats.tstat),'], paired t-test']) ;
suptitle('FIGURE 3') ; 






