clear all ; close all ; 
subs = {'sub_lyes','sub_amal','sub_lyndis','sub_valerie','sub_reihaneh','sub_samuel','sub_cesar','sub_angelina','sub_russell','sub_esteban','sub_raphael','sub_greg','sub_felix','sub_pascal'};

obliques = [40:50,130:140];
cardinals = [1:5,85:95,175:180];
verticals = [80:100];
horizontals = [1:10,170:180]; 

cd e:/saved
mdb = load('mdb'); mdb = mdb.mdb; 
mgammahigh = smooth(squeeze(mean(mean(mean(mdb(:,1:2,14:30,:),1),3),2))); 
malphalow = smooth(squeeze(mean(mean(mean(mdb(:,4:5,5:8,:),1),3),2))); 

cd e:/orientation_retinotopy/savevox ; 
bet = load_untouch_nii('bet_mni.nii.gz'); 
brainvox = find(bet.img > 0); 

for sb=1:length(subs) ; disp(subs{sb}); 

    cd(['e:\orientation_retinotopy\',subs{sb}]);    
    highc = load_untouch_nii('fnirt_mni_highcontrast_180.nii');  
    highc = reshape(highc.img,[numel(highc.img(:,:,:,1)),180]);
    highcorrs = corr(highc(brainvox,:)',mgammahigh); 
    boldh_gammacorrs = zeros(size(bet.img)); boldh_gammacorrs(brainvox) = highcorrs; 
    lowcorrs = corr(highc(brainvox,:)',malphalow); 
    boldh_alphacorrs = zeros(size(bet.img)); boldh_alphacorrs(brainvox) = lowcorrs; 
    clear highc; 
    allcorrs(:,:,:,sb,1) = boldh_gammacorrs;
    allcorrs(:,:,:,sb,2) = boldh_alphacorrs; 
    
    lowc = load_untouch_nii('fnirt_mni_lowcontrast_180.nii');     
    lowc = reshape(lowc.img,[numel(lowc.img(:,:,:,1)),180]);
    highcorrs = corr(lowc(brainvox,:)',mgammahigh); 
    boldl_gammacorrs = zeros(size(bet.img)); boldl_gammacorrs(brainvox) = highcorrs; 
    lowcorrs = corr(lowc(brainvox,:)',malphalow); 
    boldl_alphacorrs = zeros(size(bet.img)); boldl_alphacorrs(brainvox) = lowcorrs; 
    clear lowc; 
    allcorrs(:,:,:,sb,3) = boldl_gammacorrs;
    allcorrs(:,:,:,sb,4) = boldl_alphacorrs; 
    
end

cd e:/orientation_retinotopy/savevox;
mcorrs = squeeze(mean(allcorrs,4)); 
bet.img = mcorrs(:,:,:,1) ; save_untouch_nii(bet,'hgammacorrs.nii.gz'); 
bet.img = mcorrs(:,:,:,2) ; save_untouch_nii(bet,'halphacorrs.nii.gz'); 
bet.img = mcorrs(:,:,:,3) ; save_untouch_nii(bet,'lgammacorrs.nii.gz'); 
bet.img = mcorrs(:,:,:,4) ; save_untouch_nii(bet,'lalphacorrs.nii.gz'); 


hgamma = find(abs(mcorrs(:,:,:,1)) > 0.25);
lalpha = find(abs(mcorrs(:,:,:,4)) > 0.25); 
both = intersect(hgamma,lalpha); 
newbet = zeros(size(bet.img)); 
newbet(hgamma) = 3;
newbet(lalpha) = 1;
newbet(both) = 2; 
bet.img = newbet; save_untouch_nii(bet,'bothcorrs.nii.gz');






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

res_subs = reshape(allcorrs,[numel(mcorrs(:,:,:,1,1)),14,4]);
res_subs(isnan(res_subs)) = 0; 
clear meanrois 
for i=1:length(atlasinds)
   meanrois(i,:,:) = mean(res_subs(atlasinds{i},:,:),1); 
end

figure,
corrtitles = {'gamma vs BOLD (@high BOLD contrast)','alpha vs BOLD (@high BOLD contrast)','gamma vs BOLD (@low BOLD contrast)','alpha vs BOLD (@low BOLD contrast)'};
for i=1:4 ; subplot(2,2,i) ; bar(squeeze(mean(meanrois(:,:,i),2))) ; set(gca,'XTick',1:9,'XTickLabel',rnames); title(corrtitles{i}); xlabel('visual area'); ylabel('correlation(rho)') ;end










