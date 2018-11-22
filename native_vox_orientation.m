clear all ; close all ; 
subs = {'sub_lyes','sub_amal','sub_lyndis','sub_valerie','sub_reihaneh','sub_samuel','sub_cesar','sub_angelina','sub_russell','sub_esteban','sub_raphael','sub_greg','sub_felix','sub_pascal'};
subdates = {'07/06/2018','07/07/2018','01/06/2018','01/06/2018','17/05/2018','17/05/2018','11/05/2018','11/05/2018','12/04/2018','26/04/2018','26/04/2018','12/04/2018','03/05/2018','17/05/2018'};

for sb=1:length(subs) ; disp(subs{sb}); 
    tr = 0.68; 

    cd(['e:\orientation_retinotopy\',subs{sb}]);
    mask = load_untouch_nii('mean_unwarped.nii.gz'); 
    atlasinds = find(mask.img>mean(mask.img(:))*2); 
    
    lowcontrast = load_untouch_nii('blur_tproj_topup_mc_orientation_2.nii.gz');
    highcontrast = load_untouch_nii('blur_tproj_topup_mc_orientation_1.nii.gz'); 

    res_lowcontrast = reshape(lowcontrast.img,[numel(lowcontrast.img(:,:,:,1)),size(lowcontrast.img,4)]); 
    res_highcontrast = reshape(highcontrast.img,[numel(highcontrast.img(:,:,:,1)),size(highcontrast.img,4)]); 

    clear low high
    low(1,:,:) = res_lowcontrast(atlasinds,:); 
    high(1,:,:) = res_highcontrast(atlasinds,:); 
  
    startind = 10/0.68 + 6; 
    inds = round(startind:startind + (8*60)/0.68); 

    low = (low(:,:,inds));
    high = (high(:,:,inds));

    tincr = size(high,3)/8; 
    icount=1;
    clear ephigh eplow
    for i=1:tincr:size(high,3)
       ephigh(:,:,:,icount) = high(:,:,floor(i):floor(i)+floor(tincr));  
       eplow(:,:,:,icount) = low(:,:,floor(i):floor(i)+floor(tincr));  
       icount=icount+1; 
    end
    ephigh = (mean(ephigh,4)); 
    eplow = (mean(eplow,4)); 

    clear res_ephigh res_eplow
    res_ephigh(1,:,:) = (imresize(squeeze(circshift(ephigh(1,:,:),-22,3)),[size(ephigh,2),360])); 
    res_eplow(1,:,:) = (imresize(squeeze(circshift(eplow(1,:,:),-22,3)),[size(ephigh,2),360])); 
    
    res_ephigh = (res_ephigh(:,:,1:180) + res_ephigh(:,:,181:end)) / 2; 
    res_eplow = (res_eplow(:,:,1:180) + res_eplow(:,:,181:end)) / 2; 

    oratlas_high = zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3),180);
    oratlas_high = reshape(oratlas_high,[numel(oratlas_high(:,:,:,1)),180]); 
    oratlas_low = zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3),180); 
    oratlas_low = reshape(oratlas_low,[numel(oratlas_low(:,:,:,1)),180]); 
    
    oratlas_high(atlasinds,:) = res_ephigh; 
    oratlas_high = reshape(oratlas_high,[size(mask.img,1),size(mask.img,2),size(mask.img,3),180]);
    oratlas_low(atlasinds,:) = res_eplow; 
    oratlas_low = reshape(oratlas_low,[size(mask.img,1),size(mask.img,2),size(mask.img,3),180]);
  
    cat180 = load_untouch_nii('cat180.nii.gz'); 
    cat180.img = oratlas_high; save_untouch_nii(cat180,'highcontrast_180.nii.gz');
    cat180.img = oratlas_low; save_untouch_nii(cat180,'lowcontrast_180.nii.gz');
end



