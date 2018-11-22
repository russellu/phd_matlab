%{
clear all ; close all ; 
cd c:/shared/mnirefs; ls 
anats = dir('nifti*') ; 
for i=1:length(anats)
    nii = load_nii(anats(i).name) ; 
    anatvols(i,:,:,:) = nii.img ; 
    nii = load_nii(['boldcorr_',anats(i).name]) ; 
    bothcorrvols(i,:,:,:) = nii.img ; 
end
%}

clear all ; close all ; 
cd c:/shared/mnirefs2; ls 
anats = dir('sub_*') ; 
for i=1:length(anats)
    nii = load_nii(anats(i).name) ; 
    anatvols(i,:,:,:) = nii.img ; 
    nii = load_nii(['corr_',anats(i).name]) ; 
    bothcorrvols(i,:,:,:) = nii.img ; 
end

zcorrs = squeeze(mean(bothcorrvols(:,:,:,z:z+5),4)) ; 
zanats = squeeze(mean(anatvols(:,:,:,z:z+5),4)) ; 
% rotate the image so its right side up
for i=1:size(zcorrs,1) ; 
   rotcorrs(i,:,:) = fliplr(rot90(squeeze(zcorrs(i,:,:)))) ;  
   rotanats(i,:,:) = mat2gray(fliplr(rot90(squeeze(zanats(i,:,:))))) ; 
end
% want to stack the z corrs along the 1st dimension
largecorrs = zeros(size(rotcorrs,2),size(rotcorrs,1)*size(rotcorrs,3)) ; 
largeanats = zeros(size(rotcorrs,2),size(rotcorrs,1)*size(rotcorrs,3)) ; 
for i=1:size(rotcorrs,1) ; 
    largecorrs(:,((i-1)*size(rotcorrs,3)+1):(i*size(rotcorrs,3))) = squeeze(rotcorrs(i,:,:)) ; 
    largeanats(:,((i-1)*size(rotcorrs,3)+1):(i*size(rotcorrs,3))) = squeeze(rotanats(i,:,:)) ; 
end
plotoverlayIntensity2D(largeanats,largecorrs,largecorrs) ;







