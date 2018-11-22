clear all ; close all ; 
cd c:/shared/claudie ; ls 
subs=dir('*') ; subs(1:2) = [] ; 
for s=1:length(subs) ; 
    disp(subs(s).name) ; 
    cd(['c:/shared/claudie/',subs(s).name]) ; 
    roi1 = load_untouch_nii('fs_labs/3mm_lh.inferiorparietal.label.nii.gz') ; 
    roi2 = load_untouch_nii('fs_labs/3mm_lh.precentral.label.nii.gz') ;
    img1 = imdilate(roi1.img,strel(ones(3,3,3))) ; 
    img2 = imdilate(roi2.img,strel(ones(3,3,3))) ; 
    pre = load_untouch_nii('bp_anat_regf_preTMS.nii.gz') ; 
    post = load_untouch_nii('bp_anat_regf_postTMS.nii.gz') ; 
    parinds = find(img1==1) ; 
    [cpx,cpy,cpz] = ind2sub(size(img1),parinds) ; 
    motinds = find(img2==1) ; 
    [cmx,cmy,cmz] = ind2sub(size(img1),motinds) ; 
    clear pre_pts post_pts pre_mts post_mts ; 
    for i=1:length(cpx)
        pre_pts(i,:) = squeeze(pre.img(cpx(i),cpy(i),cpz(i),:)) ; 
        post_pts(i,:) = squeeze(post.img(cpx(i),cpy(i),cpz(i),:)) ; 
    end
    for i=1:length(cmx)
        pre_mts(i,:) = squeeze(pre.img(cmx(i),cmy(i),cmz(i),:)) ; 
        post_mts(i,:) = squeeze(post.img(cmx(i),cmy(i),cmz(i),:)) ; 
    end
    cpre = corr(pre_mts',pre_pts') ; 
    cpost = corr(post_mts',post_pts') ; 
    allpre{s} = cpre ; 
    allpost{s} = cpost ; 
    
    img1(parinds) = mean(cpre,1) ;
    img1(motinds) = mean(cpre,2) ;
    roi1.img = img1 ; save_untouch_nii(roi1,'roi1.nii.gz') ; 

end



%{
for i=1:11 ; 
    meanpost(i) = mean(mean(allpost{i})) ; 
    meanpre(i) = mean(mean(allpre{i})) ; 
end
meanboth = [meanpre;meanpost] ; 
[h,p,ci,tstat] = ttest(meanpre,meanpost) ; 

subplot(2,2,1) ; 
barwitherr(std(meanboth,0,2)./sqrt(11),mean(meanboth,2)) ; set(gca,'XTickLabel',{'pre','post'}) ; ylabel('mean correlation (r)') ; 
subplot(2,2,2) ; 
bar(meanboth') ; title(['p(pre>post) = ', num2str(p), ' t=',num2str(tstat.tstat)]) ;
%}





