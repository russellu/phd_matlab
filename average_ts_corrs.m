clear all ; close all ; 
subs = {'sub_lyes','sub_amal','sub_lyndis','sub_valerie','sub_reihaneh','sub_samuel','sub_cesar','sub_angelina','sub_russell','sub_esteban','sub_raphael','sub_greg','sub_felix','sub_pascal'};

cd e:/saved
mdb = load('mdb'); mdb = mdb.mdb; 
mgammahigh = smooth(squeeze(mean(mean(mean(mdb(:,1:2,14:30,:),1),3),2))); 
malphalow = smooth(squeeze(mean(mean(mean(mdb(:,4:5,5:8,:),1),3),2))); 

cd e:/ProbAtlas_v4/subj_vol_all
atlas = load_untouch_nii('maxprob_2mm.nii.gz'); 
roinames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','VO1','VO2','PHC1','PHC2','MST','hMT','LO2','LO1','V3b','V3a','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5','SPL1','FEF'};
roi_inds = {[1,2],[3,4],[5,6],[7],[8,9],[10,11],[12],[13],[14,15]};
rnames = {'V1','V2','V3','hV4','VO','PH','MST','hMT','LO'};
clear atlasinds
for i=1:length(roi_inds)
    if length(roi_inds{i}) == 2
        atlasinds{i} = [find(atlas.img == roi_inds{i}(1));find(atlas.img == roi_inds{i}(2))]; 
    else
        atlasinds{i} = find(atlas.img == roi_inds{i}(1)); 
    end
end

mni = load_untouch_nii('bet_2mm.nii.gz'); 
maskimg = imdilate(mni.img>0,strel(ones(3,3,3))); 
maskvox = find(maskimg>0); 

for sb=1:length(subs) ; disp(subs{sb}); 

    cd(['e:\orientation_retinotopy\',subs{sb}]);    
    highc = load_untouch_nii('res_fnirt_mni_highcontrast_180.nii');  
    if sb==1 ; mhigh = highc.img ; else mhigh = mhigh + highc.img; end
    lowc = load_untouch_nii('res_fnirt_mni_lowcontrast_180.nii');   
    if sb==1 ; mlow = lowc.img ; else mlow = mlow + lowc.img; end
    
    res_highc = reshape(highc.img,[numel(highc.img(:,:,:,1)),180]);
    res_lowc = reshape(lowc.img,[numel(highc.img(:,:,:,1)),180]);
    
    sb_reshighc(:,:,sb) = res_highc(maskvox,:); 
    sb_reslowc(:,:,sb) = res_lowc(maskvox,:); 
    
    for i=1:length(atlasinds)
       sub_roits(sb,i,1,:) = squeeze(mean(res_highc(atlasinds{i},:),1));  
       sub_roits(sb,i,2,:) = squeeze(mean(res_lowc(atlasinds{i},:),1));  
    end
    
end

clear corr_gammahigh corr_alphalow
[corr_gammahigh(:),p_gammahigh] = corr(mean(sb_reshighc(:,:,:),3)',mgammahigh); 
[corr_gammalow(:),p_gammalow] = corr(mean(sb_reslowc(:,:,:),3)',mgammahigh); 
[corr_alphalow(:),p_alphalow] = corr(mean(sb_reslowc(:,:,:),3)',malphalow); 
[corr_alphahigh(:),p_alphahigh] = corr(mean(sb_reshighc(:,:,:),3)',malphalow); 

cd e:/orientation_retinotopy/savevox; template = load_untouch_nii('2mm_mni.nii.gz'); 
zvals = zeros(size(highc.img(:,:,:,1))); 
zvals(maskvox) = corr_gammahigh > 0.6; ghigh = zvals; 
zvals(maskvox) = corr_gammalow > 0.6; glow = zvals; 
zvals(maskvox) = corr_alphalow< -0.6; alow = zvals; 
zvals(maskvox) = corr_alphahigh< -0.6; ahigh = zvals; 

highs = zeros(size(zvals)); highs(find(ghigh==1)) = 3; highs(find(ahigh==1)) = 1; highs(intersect(find(ahigh==1), find(ghigh==1))) = 2; 
lows = zeros(size(zvals)); lows(find(glow==1)) = 3; lows(find(alow==1)) = 1; lows(intersect(find(alow==1), find(glow==1))) = 2; 

template.img = highs ; save_untouch_nii(template,'highcorrs_2mm.nii.gz');
template.img = lows ; save_untouch_nii(template,'lowcorrs_2mm.nii.gz') ;

gsi = find(corr_gammahigh>0.65); 
asi = find(corr_alphalow<-0.65); 
bsi = intersect(gsi,asi); 

subplot(2,2,1);
plot(squeeze(mean(mean(sb_reshighc(gsi,:,:),1),3)),'Color',[1,0,0],'LineWidth',3); hold on ;
plot(squeeze(mean(mean(sb_reslowc(gsi,:,:),1),3)),'Color',[0.3,0,0],'LineWidth',3); legend('100% contrast','5% contrast'); title('gamma ROI'); xlim([0,180]); xlabel('orientation(deg)'); ylabel('BOLD (AU)'); 
subplot(2,2,2); 
plot(squeeze(mean(mean(sb_reshighc(asi,:,:),1),3)),'Color',[0,0,1],'LineWidth',3); hold on ;
plot(squeeze(mean(mean(sb_reslowc(asi,:,:),1),3)),'Color',[0,0,0.3],'LineWidth',3); legend('100% contrast','5% contrast'); title('alpha ROI'); xlim([0,180]); xlabel('orientation(deg)'); ylabel('BOLD (AU)'); 
subplot(2,2,3); 
plot(squeeze(mean(mean(sb_reshighc(bsi,:,:),1),3)),'Color',[0,1,0],'LineWidth',3); hold on ;
plot(squeeze(mean(mean(sb_reslowc(bsi,:,:),1),3)),'Color',[0,0.3,0],'LineWidth',3); legend('100% contrast','5% contrast'); title('both ROI'); xlim([0,180]); xlabel('orientation(deg)'); ylabel('BOLD (AU)'); 

newatlasimg = atlas.img(maskvox); 
clear newatlasinds
for i=1:length(roi_inds)
    if length(roi_inds{i}) == 2
        newatlasinds{i} = [find(newatlasimg == roi_inds{i}(1));find(newatlasimg == roi_inds{i}(2))]; 
    else
        newatlasinds{i} = find(newatlasimg == roi_inds{i}(1)); 
    end
end

for i=1:length(newatlasinds)
    gamma_atlasinds(i) = length(intersect(newatlasinds{i},gsi));%/length(newatlasinds{i}); 
    alpha_atlasinds(i) = length(intersect(newatlasinds{i},asi));%/length(newatlasinds{i}); 
    both_atlasinds(i) = length(intersect(newatlasinds{i},bsi));%/length(newatlasinds{i}); 
end

plot(gamma_atlasinds([1:4,9]),'r','LineWidth',3); hold on ;
plot(alpha_atlasinds([1:4,9]),'b','LineWidth',3); 
plot(both_atlasinds([1:4,9]),'g','LineWidth',3); 
set(gca,'XTick',1:5,'XTickLabel',rnames([1:4,9])); ylabel('# voxels'); xlabel('visual ROI');




for i=1:9
    [hgcorrs(i),hgps(i)] = corr(squeeze(mean(sub_roits(:,i,1,:),1)),mgammahigh); 
    [lgcorrs(i),lgps(i)] = corr(squeeze(mean(sub_roits(:,i,2,:),1)),mgammahigh); 

    [hacorrs(i),haps(i)] = corr(squeeze(mean(sub_roits(:,i,1,:),1)),malphalow); 
    [lacorrs(i),laps(i)] = corr(squeeze(mean(sub_roits(:,i,2,:),1)),malphalow); 
end

roi_colors = {[153,32,102]/255,[0,0,255]/255,[0,153,255]/255,[0,255,255]/255,[0,255,0]/255,[255,255,0]/255,[255,105,0]/255,[255,0,0]/255,[204,16,51]/255}; 

subplot(2,2,1); 
for i=1:9
   b = bar(i,hgcorrs(i)); hold on; 
   b.FaceColor = roi_colors{i}; 
end
title('gamma vs BOLD (100% contrast), *p<0.05'); ylim([-1,1]);
set(gca,'XTick',1:9,'XTickLabel',rnames,'XTickLabelRotation',45); 
for i=1:9 ; if hgps(i) < 0.05 ; text(i,.9,'*'); end; end ; ylabel('rho'); 

subplot(2,2,2); 
for i=1:9
   b = bar(i,lgcorrs(i)); hold on; 
   b.FaceColor = roi_colors{i}; 
end
title('gamma vs BOLD (5% contrast), *p<0.05'); ylim([-1,1]);
set(gca,'XTick',1:9,'XTickLabel',rnames,'XTickLabelRotation',45); 
for i=1:9 ; if lgps(i) < 0.05 ; text(i,.9,'*'); end; end ; ylabel('rho'); 

subplot(2,2,3); 
for i=1:9
   b = bar(i,hacorrs(i)); hold on; 
   b.FaceColor = roi_colors{i}; 
end
title('alpha vs BOLD (100% contrast), *p<0.05');  ylim([-1,1]);
set(gca,'XTick',1:9,'XTickLabel',rnames,'XTickLabelRotation',45); 
for i=1:9 ; if haps(i) < 0.05 ; text(i,.9,'*'); end; end; ylabel('rho'); 

subplot(2,2,4); 
for i=1:9
   b = bar(i,lacorrs(i)); hold on; 
   b.FaceColor = roi_colors{i}; 
end
title('alpha vs BOLD (5% contrast), *p<0.05');  ylim([-1,1]);
set(gca,'XTick',1:9,'XTickLabel',rnames,'XTickLabelRotation',45); 
for i=1:9 ; if laps(i) < 0.05 ; text(i,.9,'*'); end; end; ylabel('rho'); 

for i=1:9
    subplot(3,3,i)
    shadedErrorBar(1:180,squeeze(mean(sub_roits(:,i,1,:),1)),squeeze(std(sub_roits(:,i,1,:),0,1))/sqrt(14),{'Color',roi_colors{i}});
    xlim([0,180]); set(gca,'XTick',[1,45,90,135,180]); title(rnames{i});xlabel('orientation (deg)');
end


for i=1:9
    subplot(3,3,i)
    shadedErrorBar(1:180,squeeze(mean(sub_roits(:,i,2,:),1)),squeeze(std(sub_roits(:,i,2,:),0,1))/sqrt(14),{'Color',roi_colors{i}});
    xlim([0,180]); set(gca,'XTick',[1,45,90,135,180]); title(rnames{i});xlabel('orientation (deg)');
end


cd e:/saved; save('sub_roits','sub_roits'); 

