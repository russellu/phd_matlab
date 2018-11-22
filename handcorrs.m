clear all ; close all ; 
sujets = {'sujet1','sujet2','sujet3','sujet4','sujet5','sujet6','sujet8','sujet9','sujet10','sujet11','sujet12'} ; 
for s=1:length(sujets) ; 
cd(['c:/shared/claudie/',sujets{s}]) ; ls 

postnii = load_untouch_nii('bp_post.nii.gz') ; 
prenii = load_untouch_nii('bp_pre.nii.gz') ; 
graynii = load_untouch_nii('3mm_cortex.nii.gz') ; 
motornii = load_untouch_nii('3mm_motor.nii.gz') ; 

hand = find(motornii.img==1) ; 
[hx,hy,hz] = ind2sub(size(motornii.img),hand); 
for i=1:length(hx)
   hs_pre(i,:) = squeeze(prenii.img(hx(i),hy(i),hz(i),:)) ; 
   hs_post(i,:) = squeeze(postnii.img(hx(i),hy(i),hz(i),:)) ; 
end

corrpre = voxcorr(prenii.img,mean(hs_pre,1)) ; 
corrpost = voxcorr(postnii.img,mean(hs_post,1)) ; 

graynii.img = corrpre ; save_untouch_nii(graynii,'corrpre.nii.gz') ; 
graynii.img = corrpost ; save_untouch_nii(graynii,'corrpost.nii.gz') ; 


end