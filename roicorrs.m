clear all ; close all ; 
sujets = {'sujet1','sujet2','sujet3','sujet4','sujet5','sujet6','sujet8','sujet9','sujet10','sujet11','sujet12'} ; 
% get the rois
cd c:/shared/ROI/ ; 
allrois = dir('3mm_*') ; 
for r=1:length(allrois)
   roi = load_untouch_nii(allrois(r).name) ;  
   rois(r,:,:,:) = roi.img ; 
end

% correlate the rois with the 
for s=1:length(sujets) ; 
cd(['c:/shared/claudie/',sujets{s}]) ; ls 

postnii = load_untouch_nii('bp_mni_post.nii.gz') ; 
prenii = load_untouch_nii('bp_mni_pre.nii.gz') ; 
graynii = load_untouch_nii('medgm_3mm.nii.gz') ; 
%motornii = load_untouch_nii('3mm_motor.nii.gz') ; 

for r=1:size(rois,1)

hand = find(squeeze(rois(r,:,:,:))==1) ; 
[hx,hy,hz] = ind2sub(size(squeeze(rois(r,:,:,:))),hand); 
for i=1:length(hx)
   hs_pre(i,:) = squeeze(prenii.img(hx(i),hy(i),hz(i),:)) ; 
   hs_post(i,:) = squeeze(postnii.img(hx(i),hy(i),hz(i),:)) ; 
end

corrpre = voxcorr(prenii.img,mean(hs_pre,1)) ; 
corrpost = voxcorr(postnii.img,mean(hs_post,1)) ; 

graynii.img = corrpre ; save_untouch_nii(graynii,[strrep(allrois(r).name,'.nii.gz','_'),'corrpre.nii.gz']) ; 
graynii.img = corrpost ; save_untouch_nii(graynii,[strrep(allrois(r).name,'.nii.gz','_'),'corrpost.nii.gz']) ; 

end

end


