clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 
mask = load_untouch_nii('c:/shared/epireg/meanmask.nii.gz'); 
maskinds = find(mask.img==1); 
clean_fmri_names = {'bp_clean_retino_allstims_01','bp_clean_retino_allstims_02','bp_clean_retino_gamma_01','bp_clean_retino_gamma_02','bp_clean_retino_movie','bp_clean_retino_rest'};


for sb=1:length(subs)
    
    cd(['e:/fmris/badger_',subs{sb},'/atlas_fmri']);
    ls 
    for cfm = 1:length(clean_fmri_names)
        fmri = load_untouch_nii([clean_fmri_names{cfm},'.nii.gz']); 
        res_fmri = reshape(fmri.img,[numel(fmri.img(:,:,:,1)),size(fmri.img,4)]);
        globsig = squeeze(mean(res_fmri(maskinds,:),1)); 

        xcorrs = zeros(length(maskinds),41); 
        for i=1:length(maskinds)
           xcorrs(i,:) = xcorr(res_fmri(maskinds(i),20:end-20),globsig(20:end-20),20,'coeff');     
        end

        xcorrbrain = zeros(numel(mask.img(:,:,:,1)),41);
        xcorrbrain(maskinds,:) = xcorrs; 
        xcorrbrain = reshape(xcorrbrain,[size(mask.img,1),size(mask.img,2),size(mask.img,3),41]);
        allbrains(:,:,:,:,cfm,sb) = xcorrbrain;  
    end
end



allbrains(isnan(allbrains)) = 0; 
cd E:\meanepis
ref41 = load_untouch_nii('ref41.nii.gz'); 
ref41.img = squeeze(mean(mean(allbrains,5),6));
save_untouch_nii(ref41,'globalsignal_xcorr.nii.gz'); 










