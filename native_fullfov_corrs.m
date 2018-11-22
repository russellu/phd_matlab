clear all ; close all ; 
subs = {'sub_lyes','sub_amal','sub_lyndis','sub_valerie','sub_reihaneh','sub_samuel','sub_cesar','sub_angelina','sub_russell','sub_esteban','sub_raphael','sub_greg','sub_felix','sub_pascal'};


for sb=1:length(subs) ; disp(subs{sb}); 
    tr = 0.68; 

    cd(['e:\orientation_retinotopy\',subs{sb}]);
    
    fullfov = load_untouch_nii('tproj_topup_mc_fullfov.nii.gz'); 
    
    quad_design = zeros(1,round(750*tr)); 
    for i=1:2:24; quad_design((i-1)*10+1:(i-1)*10+10) = 1; end
    for i=1:2:24; quad_design((i-1)*10+1+250:(i-1)*10+10+250) = 1; end
    quad_design = circshift(quad_design,6); quad_design = smooth(imresize(quad_design,[1,750])); 
    fullfov_design = quad_design; 
    h1_fullfov = fullfov.img(:,:,:,1:end/2); 
    h2_fullfov = fullfov.img(:,:,:,376:end); 
    h1_fullfov_design = fullfov_design(1:end/2);
    h2_fullfov_design = fullfov_design(376:end); 
    
    corrs1_fullfov = voxcorr(h1_fullfov(:,:,:,20:end-20),h1_fullfov_design(20:end-20)); 
    corrs2_fullfov = voxcorr(h2_fullfov(:,:,:,20:end-20),h2_fullfov_design(20:end-20)); 
    
    template = load_untouch_nii('mean_unwarped.nii.gz');
    template.img = corrs2_fullfov; 
    save_untouch_nii(template,'fovcorrs.nii.gz'); 
       
end






