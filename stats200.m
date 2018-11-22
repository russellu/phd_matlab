clear all ; close all ; 
sujets = {'sujet1','sujet2','sujet3','sujet4','sujet5','sujet6','sujet8','sujet9','sujet10','sujet11','sujet12'} ; 

for s=1:length(sujets) ; 
cd(['c:/shared/claudie/',sujets{s}]) ; ls 
pre = load_untouch_nii('corrpre.nii.gz') ;
post = load_untouch_nii('corrpost.nii.gz') ; 
allpre(s,:,:,:,:) = pre.img ; 
allpost(s,:,:,:,:) = post.img ; 
end


nii = load_untouch_nii('mnicat1.nii.gz') ;
for roi=1:200 ; disp(roi) ; 
vecpre = reshape(squeeze(allpre(:,:,:,:,roi)),[11,numel(allpre(1,:,:,:,roi))]) ; 
vecpost = reshape(squeeze(allpost(:,:,:,:,roi)),[11,numel(allpost(1,:,:,:,roi))]) ; 
vecpost(isnan(vecpost)) = 0 ; vecpre(isnan(vecpre)) = 0 ; vecpost(isinf(vecpost)) = 0 ; vecpost(isinf(vecpost)) = 0 ; 
[t,p,ci,stats] = ttest(vecpost,vecpre) ; 
tbrain = reshape(stats.tstat,size(nii.img)) ; 
max(max(max(tbrain)))
pbrain = reshape(p,size(nii.img)) ; 
nii.img = (tbrain) ; save_untouch_nii(nii,['finalts/t_',num2str(roi),'.nii.gz']) 
nii.img = (pbrain) ; save_untouch_nii(nii,['finalts/p_',num2str(roi),'.nii.gz']) 
end