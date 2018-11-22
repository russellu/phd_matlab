clear all ; close all ; 
sujets = {'sujet1','sujet2','sujet3','sujet4','sujet5','sujet6','sujet8','sujet9','sujet10','sujet11','sujet12'} ; 

for s=1:length(sujets) ; 
cd(['c:/shared/claudie/',sujets{s}]) ; ls 
nii = load_untouch_nii('den_alff_pre_ALFF.nii.gz') ; 
allpres(s,:,:,:) = nii.img ; 
nii = load_untouch_nii('den_alff_post_ALFF.nii.gz') ; 
allpost(s,:,:,:) = nii.img ; 
end

vecpre = reshape(allpres,[11,numel(allpres(1,:,:,:))]) ; 
vecpost = reshape(allpost,[11,numel(allpost(1,:,:,:))]) ; 
vecpost(isnan(vecpost)) = 0 ; vecpre(isnan(vecpre)) = 0 ; vecpost(isinf(vecpost)) = 0 ; vecpost(isinf(vecpost)) = 0 ; 
[t,p,ci,stats] = ttest(vecpost,vecpre) ; 
tbrain = reshape(stats.tstat,size(nii.img)) ; 
pbrain = reshape(p,size(nii.img)) ; 
nii.img = double(tbrain) ; save_untouch_nii(nii,'t_denoised_ALFF.nii.gz') 
nii.img = double(pbrain) ; save_untouch_nii(nii,'p_denoised_ALFF.nii.gz') 

dilt = imdilate(tbrain>6,strel(ones(9,9,9))) ; 
bw = bwconncomp(dilt) ; 
sz = cellfun(@length,bw.PixelIdxList) ; 
maxclust = find(sz==max(sz)) ; 
inds = bw.PixelIdxList{maxclust} ; 
zs = zeros(size(dilt)) ; zs(inds) = 1 ; zs = imerode(zs,strel(ones(5,5,5))) ; 
for i=1:11
    prebrain = squeeze(allpres(i,:,:,:)) ; 
    postbrain = squeeze(allpost(i,:,:,:)) ; 
    mpres(i) = squeeze(mean(prebrain(zs==1))) ; 
    mpost(i) = squeeze(mean(postbrain(zs==1))) ; 
end

both = [mpres;mpost] ; 
bar(both') ; xlabel('subject #') ; ylabel('reho') ; legend({'pre TMS','post TMS'}) ; 
[h,p,ci,stats]= ttest(both(2,:),both(1,:)) ; 

barwitherr(squeeze(std(both,0,2))./sqrt(11),mean(both,2)) ; set(gca,'XTickLabel',{'pre','post'}) ; ylabel('reho') ; 
title(['paired t=',num2str(stats.tstat),', p=',num2str(p)]) ; 







