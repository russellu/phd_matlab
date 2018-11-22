clear all ;close all ; 
subs = {'a_alex','a_charest','a_esteban','a_fabio','a_gab','a_gabriella','a_genevieve','a_gina','a_guillaume','a_jeremie','a_julie','a_katrine','a_lisa','a_marc',...
    'a_marie','a_mathieu','a_maxime','a_mingham','a_patricia','a_po','a_russell','a_sunachakan','a_tah','a_vincent'} ; 

for sb=1:length(subs)
    
 
    cd(['e:/nimg_pool/',subs{sb}]) ; 
    nii = load_untouch_nii('alphabeta_mni.nii.gz'); 
    nii.img(isnan(nii.img)) = 0; 
    alpha_nii(:,:,:,sb) = nii.img; 
    
    nii = load_untouch_nii('gamma_mni.nii.gz'); 
    nii.img(isnan(nii.img)) = 0; 
    gamma_nii(:,:,:,sb) = nii.img; 
    
    corrs = load_untouch_nii('cleancorrs_mni.nii.gz'); 
    corr_nii(:,:,:,sb) = corrs.img; 
end

cd c:/shared/ATLASES ;
mni = load_untouch_nii('MNI152_T1_1mm.nii.gz');



gamma_olap = zeros(size(nii.img)); 
alphabeta_olap = zeros(size(nii.img)); 
corr_olap = zeros(size(nii.img)); 

for i=1:24
   sub_alpha = alpha_nii(:,:,:,i); 
   [sv,si] = sort(sub_alpha(:),'ascend');
   alphabeta_olap(si(1:15000)) = alphabeta_olap(si(1:15000)) + 1; 
   
   sub_gamma = gamma_nii(:,:,:,i); 
   [sv,si] = sort(sub_gamma(:),'descend');
   gamma_olap(si(1:15000)) = gamma_olap(si(1:15000)) + 1; 
   
   sub_corr = corr_nii(:,:,:,i); 
   [sv,si] = sort(sub_corr(:),'descend');
   corr_olap(si(1:15000)) = corr_olap(si(1:15000)) + 1; 
    
end

alphabeta_olap = alphabeta_olap/24; 
gamma_olap = gamma_olap/24;
corr_olap = corr_olap/24; 

nii.img = alphabeta_olap; save_untouch_nii(nii,'sf_alphabeta_olap.nii.gz'); 
nii.img = gamma_olap; save_untouch_nii(nii,'sf_gamma_olap.nii.gz'); 
nii.img = corr_olap; save_untouch_nii(nii,'sf_corr_olap.nii.gz'); 

nii.img = squeeze(mean(gamma_nii,4))/1000000000;  save_untouch_nii(nii,'sf_gamma_amp.nii.gz'); 
nii.img = squeeze(mean(alpha_nii,4))/1000000000;  save_untouch_nii(nii,'sf_alphabeta_amp.nii.gz'); 
nii.img = squeeze(mean(corr_nii,4));  save_untouch_nii(nii,'sf_corr_amp.nii.gz'); 








