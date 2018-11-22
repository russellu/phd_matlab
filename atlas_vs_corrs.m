clear all ; close all ; 
cd E:\meanepis

%atlas = load_untouch_nii('atlas_vein.nii.gz'); 
%res_atlas = atlas.img(:); 

atlas = load_untouch_nii('motion_xcorrs.nii.gz'); 
res_atlas = atlas.img(:,:,:,28); 
res_atlas = res_atlas(:); 
meanmask = load_untouch_nii('c:/shared/epireg/meanmask.nii.gz'); 

corrs = load_untouch_nii('mean_all_vox.nii.gz'); 
res_corrimg = reshape(corrs.img,[numel(corrs.img(:,:,:,1)),size(corrs.img,4)]); 

maskinds = find(meanmask.img==1); 


clear rhos ps
for i=1:50
    vals1 = res_corrimg(maskinds,i); 
    vals2 = res_atlas(maskinds); 
    [rhos(i),ps(i)] = corr(res_corrimg(maskinds(vals1~=0 & vals2~=0),i),res_atlas(maskinds(vals1~=0 & vals2 ~= 0)));
end


atlasimg = atlas.img(:,:,:,28); 
subplot(1,2,1) ; imagesc(imrotate(squeeze(mean(atlasimg(:,:,8:15),3)),270)); colormap jet; 
corrimg = mean(corrs.img(:,:,:,30:end),4); 
subplot(1,2,2) ; imagesc(imrotate(squeeze(mean(corrimg(:,:,8:15),3)),270)); colormap jet; 




