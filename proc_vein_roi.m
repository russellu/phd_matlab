cd E:\rawbadger\badger_mri\alex\nii

vein = load_untouch_nii('vein_in_epi.nii.gz'); 
fmri = load_untouch_nii('bp_reg_topup_mc_retino_rest.nii.gz'); 
vein_inds = find(vein.img>0); 
res_fmri = reshape(fmri.img,[numel(fmri.img(:,:,:,1)),size(fmri.img,4)]); 

corrs = voxcorr(fmri.img(:,:,:,20:end-20),squeeze(mean(res_fmri(vein_inds,20:end-20),1))); 

for i=1:33
   subplot(3,11,i) ; 
   imagesc(imrotate(squeeze(corrs(:,:,i)),270),[-.5,.5]) ; colormap jet; 
end

fref = load_untouch_nii('fref.nii.gz'); 
fref.img = corrs; 
save_untouch_nii(fref,'veincorrs.nii.gz'); 


vein_vox = squeeze(mean(res_fmri(vein_inds,:),1)); 

fmri_xcorrs = zeros(size(res_fmri,1),41); 
for i=1:size(fmri_xcorrs,1)
   fmri_xcorrs(i,:) = xcorr(res_fmri(i,20:end-20),vein_vox(20:end-20),20,'coeff'); 
    
end

fmri_xcorrs = reshape(fmri_xcorrs,[size(fmri.img,1),size(fmri.img,2),size(fmri.img,3),41]); 
save_nii(make_nii(fmri_xcorrs),'vein_xcorrs.nii.gz'); 

