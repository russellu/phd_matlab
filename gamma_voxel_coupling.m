clear all; close all; 
cd E:\sim\saved
%  freqzimgs(:,:,:,fr,sc,sb) = zimg; 

corrs = load('mc_yes_gsr'); corrs = corrs.mc_yes_gsr; 
temp = zeros(size(corrs,1),size(corrs,2),size(corrs,3),size(corrs,4),4,size(corrs,6)); 
temp(:,:,:,:,1,:) = squeeze(mean(corrs(:,:,:,:,1:2,:),5)); 
temp(:,:,:,:,2,:) = squeeze(mean(corrs(:,:,:,:,3:4,:),5)); 
temp(:,:,:,:,3,:) = squeeze(mean(corrs(:,:,:,:,5,:),5)); 
temp(:,:,:,:,4,:) = squeeze(mean(corrs(:,:,:,:,6,:),5)); 
corrs = temp; 


corrs2 = load('mc_no_gsr'); corrs2 = corrs2.mc_no_gsr; 
temp = zeros(size(corrs2,1),size(corrs2,2),size(corrs2,3),size(corrs2,4),4,size(corrs2,6)); 
temp(:,:,:,:,1,:) = squeeze(mean(corrs2(:,:,:,:,1:2,:),5)); 
temp(:,:,:,:,2,:) = squeeze(mean(corrs2(:,:,:,:,3:4,:),5)); 
temp(:,:,:,:,3,:) = squeeze(mean(corrs2(:,:,:,:,5,:),5)); 
temp(:,:,:,:,4,:) = squeeze(mean(corrs2(:,:,:,:,6,:),5)); 
corrs2 = temp; 


slices = [5,8,10,12,14,16,20]; 
cd E:\sim\fmri\alex ; meanatlas = load_untouch_nii('mean.nii.gz'); 

thresh = 0.1; 
for sl=1:length(slices)
    subplottight(4,7,sl); 
    corrimg = squeeze(mean(mean(corrs(:,:,slices(sl),5,2,:),3),6)); corrimg(1,1) = -thresh; corrimg(1,2) = thresh; corrimg(corrimg>thresh) = thresh; corrimg(corrimg<-thresh) = -thresh; 
    corrimg(isnan(corrimg)) = 0; 
    plotoverlayIntensity2D(meanatlas.img(:,:,slices(sl)),mat2gray(abs(corrimg)),corrimg,270);    
end

thresh = 0.1; 
for sl=1:length(slices)
    subplottight(4,7,sl+7); 
    corrimg = squeeze(mean(mean(corrs2(:,:,slices(sl),5,2,:),3),6)); corrimg(1,1) = -thresh; corrimg(1,2) = thresh; corrimg(corrimg>thresh) = thresh; corrimg(corrimg<-thresh) = -thresh; 
    corrimg(isnan(corrimg)) = 0; 
    plotoverlayIntensity2D(meanatlas.img(:,:,slices(sl)),mat2gray(abs(corrimg)),corrimg,270);    
end

thresh = 0.1; 
for sl=1:length(slices)
    subplottight(4,7,sl+14); 
    corrimg = squeeze(mean(mean(corrs(:,:,slices(sl),5,4,:),3),6)); corrimg(1,1) = -thresh; corrimg(1,2) = thresh; corrimg(corrimg>thresh) = thresh; corrimg(corrimg<-thresh) = -thresh; 
    corrimg(isnan(corrimg)) = 0; 
    plotoverlayIntensity2D(meanatlas.img(:,:,slices(sl)),mat2gray(abs(corrimg)),corrimg,270);    
end

thresh = 0.1; 
for sl=1:length(slices)
    subplottight(4,7,sl+21); 
    corrimg = squeeze(mean(mean(corrs2(:,:,slices(sl),5,4,:),3),6)); corrimg(1,1) = -thresh; corrimg(1,2) = thresh; corrimg(corrimg>thresh) = thresh; corrimg(corrimg<-thresh) = -thresh; 
    corrimg(isnan(corrimg)) = 0; 
    plotoverlayIntensity2D(meanatlas.img(:,:,slices(sl)),mat2gray(abs(corrimg)),corrimg,270);    
end






