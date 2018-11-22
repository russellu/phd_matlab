clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','sukhman','valerie'}; 
scans = {'bp_reg_topup_mc_retino_allstims_01','bp_reg_topup_mc_retino_allstims_02','bp_reg_topup_mc_retino_gamma_01','bp_reg_topup_mc_retino_gamma_02','bp_reg_topup_mc_retino_movie','bp_reg_topup_mc_retino_rest'}; 

for sb=1:length(subs)
    cd(['E:\rawbadger\badger_mri\',subs{sb},'\nii']);   
    
    seed = load_untouch_nii('seed_epi.nii.gz'); 
    binseed = seed.img>0.1; 
    clusts = bwconncomp(binseed); 
    clusts.NumObjects;
    
    clustmask = zeros(size(binseed)); 
    clustmask(clusts.PixelIdxList{1}) = 1; 
    clustmask(clusts.PixelIdxList{2}) = 2; 
    clustmask = imdilate(clustmask,strel(ones(3,3,3))); 
    
    f41= load_untouch_nii('f_41.nii.gz'); 
%figure
    for sc=1:length(scans)
        fmri = load_untouch_nii(['warped/',scans{sc},'.nii.gz']); 
        res_fmri = reshape(fmri.img,[numel(fmri.img(:,:,:,1)),size(fmri.img,4)]); 
        m1 = zscore(mean(res_fmri(find(clustmask==1),:),1)); 
        m2 = zscore(mean(res_fmri(find(clustmask==2),:),1)); 
 %       subplot(2,3,sc); plot(m1) ; hold on ;plot(m2)
        xcorrs(sb,sc,:) = xcorr((m2(20:end-20)),(m1(20:end-20)),20,'coeff'); 
        
        fcorrs_1 = zeros(size(res_fmri,1),41); 
        for i=1:size(fcorrs_1,1)
           fcorrs_1(i,:) = xcorr(res_fmri(i,20:end-20),m1(20:end-20),20,'coeff');  
        end
        
        res_fcorrs_1 = reshape(fcorrs_1,[size(fmri.img,1),size(fmri.img,2),size(fmri.img,3),41]); 
        f41.img = res_fcorrs_1;
        save_untouch_nii(f41,['warp_fcorr2_',scans{sc},'.nii.gz']);
        
        fcorrs_2 = zeros(size(res_fmri,1),41); 
        for i=1:size(fcorrs_2,1)
           fcorrs_2(i,:) = xcorr(res_fmri(i,20:end-20),m2(20:end-20),20,'coeff');  
        end
        
        res_fcorrs_2 = reshape(fcorrs_2,[size(fmri.img,1),size(fmri.img,2),size(fmri.img,3),41]); 
        f41.img = res_fcorrs_2;
        save_untouch_nii(f41,['warp_fcorr1_',scans{sc},'.nii.gz']);
        
    end
     
end
%imagesc(squeeze(mean(xcorrs,2)))
