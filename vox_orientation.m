clear all ; close all ; 
subs = {'sub_lyes','sub_amal','sub_lyndis','sub_valerie','sub_reihaneh','sub_samuel','sub_cesar','sub_angelina','sub_russell','sub_esteban','sub_raphael','sub_greg','sub_felix','sub_pascal'};
subdates = {'07/06/2018','07/07/2018','01/06/2018','01/06/2018','17/05/2018','17/05/2018','11/05/2018','11/05/2018','12/04/2018','26/04/2018','26/04/2018','12/04/2018','03/05/2018','17/05/2018'};

cd E:\orientation_retinotopy\fatlas
atlas = load_untouch_nii('mean.nii.gz'); 
binatlas = double(atlas.img > mean(atlas.img(:))*2); 
atlasinds = find(binatlas==1); 

v1atlas = load_untouch_nii('v1atlas.nii.gz'); 

for sb=1:length(subs) ; disp(subs{sb}); 
    tr = 0.68; 

    cd(['e:\orientation_retinotopy\',subs{sb}]);
    
    lowcontrast = load_untouch_nii('reg_bp_mc_orientation_2.nii');
    highcontrast = load_untouch_nii('reg_bp_mc_orientation_1.nii'); 

    res_lowcontrast = reshape(lowcontrast.img,[numel(lowcontrast.img(:,:,:,1)),size(lowcontrast.img,4)]); 
    res_highcontrast = reshape(highcontrast.img,[numel(highcontrast.img(:,:,:,1)),size(highcontrast.img,4)]); 
    
    low(sb,:,:) = res_lowcontrast(atlasinds,:); 
    high(sb,:,:) = res_highcontrast(atlasinds,:); 
  
end

startind = 10/0.68 + 6; 
inds = round(startind:startind + (8*60)/0.68); 

low = squeeze(low(:,:,inds));
high = squeeze(high(:,:,inds));

tincr = size(high,3)/8; 
icount=1;
clear ephigh eplow
for i=1:tincr:size(high,3)
   ephigh(:,:,:,icount) = high(:,:,floor(i):floor(i)+floor(tincr));  
   eplow(:,:,:,icount) = low(:,:,floor(i):floor(i)+floor(tincr));  
   icount=icount+1; 
end
ephigh = squeeze(mean(ephigh,4)); 
eplow = squeeze(mean(eplow,4)); 

for i=1:14
   res_ephigh(i,:,:) = (imresize(squeeze(circshift(ephigh(i,:,:),-22,3)),[size(ephigh,2),360])); 
   res_eplow(i,:,:) = (imresize(squeeze(circshift(eplow(i,:,:),-22,3)),[size(ephigh,2),360])); 
end

res_ephigh = (res_ephigh(:,:,1:180) + res_ephigh(:,:,181:end)) / 2; 
res_eplow = (res_eplow(:,:,1:180) + res_eplow(:,:,181:end)) / 2; 

cd e:/saved
mdb = load('mdb'); mdb = mdb.mdb; 
mgammahigh = smooth(squeeze(mean(mean(mean(mdb(:,1:2,14:30,:),1),3),2))); 
malphalow = smooth(squeeze(mean(mean(mean(mdb(:,4:5,5:8,:),1),3),2))); 
m_ephigh = squeeze(mean(res_ephigh,1)); 
m_eplow = squeeze(mean(res_eplow,1)); 
for i=1:size(m_ephigh,1)
   highcorrs_gamma(i) = corr2(mgammahigh',m_ephigh(i,:));  
   lowcorrs_gamma(i) = corr2(mgammahigh',m_eplow(i,:));    
   highcorrs_alpha(i) = corr2(malphalow',m_ephigh(i,:));  
   lowcorrs_alpha(i) = corr2(malphalow',m_eplow(i,:)); 
end

gamma_atlas = zeros(size(binatlas));
alpha_atlas = zeros(size(binatlas)); 

zatlas = zeros(size(binatlas)); 
[sv,si] = sort(highcorrs_gamma.^2,'descend'); 
zatlas(atlasinds(si(1:1000))) = 3; 
gamma_atlas(atlasinds) = highcorrs_gamma; 

[sv,si] = sort(lowcorrs_alpha.^2,'descend'); 
zatlas(atlasinds(si(1:1000))) = zatlas(atlasinds(si(1:1000))) - 1 ; 
alpha_atlas(atlasinds) = lowcorrs_alpha; 

newzatlas = zeros(size(zatlas)); 
newzatlas(zatlas==-1) = 1; 
newzatlas(zatlas==2) = 2;
newzatlas(zatlas==3) = 3; 

atlas.img = newzatlas; 
save_untouch_nii(atlas,'gamma_alpha_topinds.nii.gz');

v1_inds = find(v1atlas.img==2);
v2_inds = find(v1atlas.img==1); 
gamma_pref = gamma_atlas.^2 > alpha_atlas.^2; 

v1_prop_gamma = sum(gamma_pref(v1_inds))/length(v1_inds); 
v2_prop_gamma = sum(gamma_pref(v2_inds))/length(v2_inds); 
v1_prop_alpha = 1-v1_prop_gamma;
v2_prop_alpha = 1-v2_prop_gamma; 

h=bar(1,v1_prop_gamma); h.FaceColor = [1,0,0]; hold on; h=bar(2,v1_prop_alpha); h.FaceColor = [0,0,1]; 
h=bar(4,v2_prop_gamma); h.FaceColor = [1,0,0]; h=bar(5,v2_prop_alpha); h.FaceColor = [0,0,1]; 









