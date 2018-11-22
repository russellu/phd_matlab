clear all ; close all ; 
subs = {'sub_lyes','sub_amal','sub_lyndis','sub_valerie','sub_reihaneh','sub_samuel','sub_cesar','sub_angelina','sub_russell','sub_esteban','sub_raphael','sub_greg','sub_felix','sub_pascal'};

obliques = [40:50,130:140];
cardinals = [1:5,85:95,175:180];
verticals = [80:100];
horizontals = [1:10,170:180]; 

for sb=1:length(subs) ; disp(subs{sb}); 

    cd(['e:\orientation_retinotopy\',subs{sb}]);
    
    highc = load_untouch_nii('fnirt_mni_highcontrast_180.nii'); 
    
    highsubimg(:,:,:,sb,1) = mean(highc.img(:,:,:,obliques),4); 
    highsubimg(:,:,:,sb,2) = mean(highc.img(:,:,:,cardinals),4); 
    highsubimg(:,:,:,sb,3) = mean(highc.img(:,:,:,verticals),4);
    highsubimg(:,:,:,sb,4) = mean(highc.img(:,:,:,horizontals),4); 
    clear highc 
    lowc = load_untouch_nii('fnirt_mni_lowcontrast_180.nii'); 
    lowsubimg(:,:,:,sb,1) = mean(lowc.img(:,:,:,obliques),4); 
    lowsubimg(:,:,:,sb,2) = mean(lowc.img(:,:,:,cardinals),4); 
    lowsubimg(:,:,:,sb,3) = mean(lowc.img(:,:,:,verticals),4);
    lowsubimg(:,:,:,sb,4) = mean(lowc.img(:,:,:,horizontals),4);    
    clear lowc
    
end

reshighsubs = reshape(highsubimg,[numel(highsubimg(:,:,:,1,1)),14,4]);
reslowsubs = reshape(lowsubimg,[numel(lowsubimg(:,:,:,1,1)),14,4]);
cd e:/orientation_retinotopy/savevox ; mni = load_untouch_nii('MNI152_T1_1mm.nii'); 

[h,p,ci,stats] = ttest(reshighsubs(:,:,1)',reshighsubs(:,:,2)'); 
res_obliquetvals = reshape(stats.tstat,size(highsubimg(:,:,:,1,1))); mni.img = res_obliquetvals; save_untouch_nii(mni,'high_oblique.nii.gz'); 
ob_pvals = reshape(p,size(highsubimg(:,:,:,1,1))); bin_ob = double(ob_pvals<0.01); 

[h,p,ci,stats] = ttest(reshighsubs(:,:,3)',reshighsubs(:,:,4)'); 
res_verticaltvals = reshape(stats.tstat,size(highsubimg(:,:,:,1,1)));mni.img = res_verticaltvals; save_untouch_nii(mni,'high_vertical.nii.gz'); 
vert_pvals = reshape(p,size(highsubimg(:,:,:,1,1))); bin_vert  =double(vert_pvals<0.01); 

bin_both = zeros(size(vert_pvals)); bin_both(bin_vert==1) = 1; bin_both(bin_ob==1) =3; bin_both(unique(find(bin_vert==1 & bin_ob==1))) = 2; 
mni.img = bin_both ; save_untouch_nii(mni,'highcont_binboth.nii.gz'); 

[h,p,ci,stats] = ttest(reslowsubs(:,:,1)',reslowsubs(:,:,2)'); 
res_obliquetvals = reshape(stats.tstat,size(highsubimg(:,:,:,1,1)));mni.img = res_obliquetvals; save_untouch_nii(mni,'low_oblique.nii.gz'); 
ob_pvals = reshape(p,size(highsubimg(:,:,:,1,1))); bin_ob = double(ob_pvals<0.01); 

[h,p,ci,stats] = ttest(reslowsubs(:,:,3)',reslowsubs(:,:,4)'); 
res_verticaltvals = reshape(stats.tstat,size(highsubimg(:,:,:,1,1)));mni.img = res_verticaltvals; save_untouch_nii(mni,'low_cardinal.nii.gz'); 
vert_pvals = reshape(p,size(highsubimg(:,:,:,1,1))); bin_vert  =double(vert_pvals<0.01); 

bin_both = zeros(size(vert_pvals)); bin_both(bin_vert==1) = 1; bin_both(bin_ob==1) =3; bin_both(unique(find(bin_vert==1 & bin_ob==1))) = 2; 
mni.img = bin_both ; save_untouch_nii(mni,'lowcont_binboth.nii.gz'); 


cd e:/ProbAtlas_v4/subj_vol_all
leftatlas = load_untouch_nii('maxprob_vol_lh.nii.gz'); 
rightatlas = load_untouch_nii('maxprob_vol_rh.nii.gz'); 
roinames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','VO1','VO2','PHC1','PHC2','MST','hMT','LO2','LO1','V3b','V3a','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5','SPL1','FEF'};
roi_inds = {[1,2],[3,4],[5,6],[7],[8,9],[10,11],[12],[13],[14,15]};
rnames = {'V1','V2','V3','hV4','VO','PH','MST','hMT','LO'};

clear atlasinds
for i=1:length(roi_inds)
    if length(roi_inds{i}) == 2
        atlasinds{i} = [find(leftatlas.img + rightatlas.img == roi_inds{i}(1));find(leftatlas.img + rightatlas.img == roi_inds{i}(2))]; 
    else
        atlasinds{i} = find(leftatlas.img + rightatlas.img == roi_inds{i}(1)); 
    end
end
clear meanhighrois meanlowrois
for i=1:length(atlasinds)
   meanhighrois(i,:,:) = mean(reshighsubs(atlasinds{i},:,:),1); 
   meanlowrois(i,:,:) = mean(reslowsubs(atlasinds{i},:,:),1); 
end
clear obhighroits obhighroips oblowroits oblowroips verthighroits verthighroips vertlowroits vertlowroips 
for i=1:length(atlasinds)
   [h,p,ci,stats] = ttest(squeeze(meanhighrois(i,:,1)),squeeze(meanhighrois(i,:,2)));  
   obhighroits(i) = stats.tstat; obhighroips(i) = p; 
   
   [h,p,ci,stats] = ttest(squeeze(meanhighrois(i,:,3)),squeeze(meanhighrois(i,:,4)));  
   verthighroits(i) = stats.tstat; verthighroips(i) = p; 
   
   [h,p,ci,stats] = ttest(squeeze(meanlowrois(i,:,1)),squeeze(meanlowrois(i,:,2)));  
   oblowroits(i) = stats.tstat; oblowroips(i) = p; 
   
   [h,p,ci,stats] = ttest(squeeze(meanlowrois(i,:,3)),squeeze(meanlowrois(i,:,4)));  
   vertlowroits(i) = stats.tstat; vertlowroips(i) = p; 
   
end


roi_colors = {[102,49,153]/255,[0,0,255]/255,[0,153,255]/255,[0,255,255]/255,[0,255,0]/255,[255,255,0]/255,[255,153,0]/255,[255,68,0]/255,[204,16,51]/255}; 


figure,
subplot(2,2,1); 
for i=1:9
   b = bar(i,obhighroits(i)); hold on; 
   b.FaceColor = roi_colors{i}; 
end
title('100% contrast, t-value, oblique - cardinal, *p<0.05'); ylim([-3,7.5]);
set(gca,'XTick',1:9,'XTickLabel',rnames,'XTickLabelRotation',45); 
for i=1:9 ; if obhighroips(i) < 0.05 ; text(i,7,'*'); end; end ; ylabel('t-value'); 

subplot(2,2,2); 
for i=1:9
   b = bar(i,verthighroits(i)); hold on; 
   b.FaceColor = roi_colors{i}; 
end
title('100% contrast, t-value, vertical - horizontal, *p<0.05'); ylim([-3,7.5]);
set(gca,'XTick',1:9,'XTickLabel',rnames,'XTickLabelRotation',45); 
for i=1:9 ; if verthighroips(i) < 0.05 ; text(i,7,'*'); end; end ; ylabel('t-value'); 

subplot(2,2,3); 
for i=1:9
   b = bar(i,oblowroits(i)); hold on; 
   b.FaceColor = roi_colors{i}; 
end
title('5% contrast, t-value, oblique - cardinal, *p<0.05');  ylim([-3,7.5]);
set(gca,'XTick',1:9,'XTickLabel',rnames,'XTickLabelRotation',45); 
for i=1:9 ; if oblowroips(i) < 0.05 ; text(i,7,'*'); end; end; ylabel('t-value'); 

subplot(2,2,4); 
for i=1:9
   b = bar(i,vertlowroits(i)); hold on; 
   b.FaceColor = roi_colors{i}; 
end
title('5% contrast, t-value, vertical - horizontal, *p<0.05');  ylim([-3,7.5]);
set(gca,'XTick',1:9,'XTickLabel',rnames,'XTickLabelRotation',45); 
for i=1:9 ; if vertlowroips(i) < 0.05 ; text(i,7,'*'); end; end; ylabel('t-value'); 





    