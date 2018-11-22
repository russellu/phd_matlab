clear all ; close all ; 
subs = {'sub_lyes','sub_amal','sub_lyndis','sub_valerie','sub_reihaneh','sub_samuel','sub_cesar','sub_angelina','sub_russell','sub_esteban','sub_raphael','sub_greg','sub_felix','sub_pascal'};
subdates = {'07/06/2018','07/07/2018','01/06/2018','01/06/2018','17/05/2018','17/05/2018','11/05/2018','11/05/2018','12/04/2018','26/04/2018','26/04/2018','12/04/2018','03/05/2018','17/05/2018'};


for sb=3%:length(subs) ; disp(subs{sb}); 
    tr = 0.68; 

    cd(['e:\orientation_retinotopy\',subs{sb}]);
    
    fullfov = load_untouch_nii('bp_fullfov.nii.gz'); 
    lowcontrast = load_untouch_nii('topup_mc_orientation_2.nii.gz');
    highcontrast = load_untouch_nii('topup_mc_orientation_1.nii.gz'); 
    
    quad_design = zeros(1,round(750*tr)); 
    for i=1:2:24; quad_design((i-1)*10+1:(i-1)*10+10) = 1; end
    for i=1:2:24; quad_design((i-1)*10+1+250:(i-1)*10+10+250) = 1; end
    quad_design = circshift(quad_design,6); quad_design = smooth(imresize(quad_design,[1,750])); 
    hemis_design = quad_design; 
    fullfov_design = quad_design; 
    h1_fullfov = fullfov.img(:,:,:,1:end/2); 
    h2_fullfov = fullfov.img(:,:,:,376:end); 
    h1_fullfov_design = fullfov_design(1:end/2);
    h2_fullfov_design = fullfov_design(376:end); 
    
    corrs1_fullfov = voxcorr(h1_fullfov(:,:,:,20:end-20),h1_fullfov_design(20:end-20)); 
    corrs2_fullfov = voxcorr(h2_fullfov(:,:,:,20:end-20),h2_fullfov_design(20:end-20)); 
    
        [sv,si] = sort(corrs1_fullfov(:),'ascend'); 

        
             
    res_lowcontrast = reshape(lowcontrast.img,[numel(lowcontrast.img(:,:,:,1)),size(lowcontrast.img,4)]); 
    res_highcontrast = reshape(highcontrast.img,[numel(highcontrast.img(:,:,:,1)),size(highcontrast.img,4)]); 
        res_fullfov = reshape(fullfov.img,[numel(fullfov.img(:,:,:,1)),size(fullfov.img,4)]); 

    
    corrvox = res_highcontrast(si(1:2000),:); 
    corrvox = eegfiltfft(corrvox,1/0.693,0.005,1); 
    [weights,sphere] = runica(corrvox','maxsteps',1024);    
    acts = weights*sphere*corrvox; 
    winv = pinv(weights*sphere); 
    
    acts = fastICA(corrvox,10,'Kurtosis'); 
    
    [u,s,v] = svd(corrvo'); 
    
    
end

