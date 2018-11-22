clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 
source_folders = {'den_retino_gamma_01','den_retino_gamma_02','den_retino_movie','den_retino_rest'}; %'den_retino_allstims_01','den_retino_allstims_02',
clean_fmri_names = {'bp_clean_retino_gamma_01','bp_clean_retino_gamma_02','bp_clean_retino_movie','bp_clean_retino_rest'}; %'bp_clean_retino_allstims_01','bp_clean_retino_allstims_02',
gsr_fmri_names = {'bp_gsr_clean_retino_gamma_01','bp_gsr_clean_retino_gamma_02','bp_gsr_clean_retino_movie','bp_gsr_clean_retino_rest'};%'bp_gsr_clean_retino_allstims_01','bp_gsr_clean_retino_allstims_02',
eeg_set_names = {'retino_gamma_01','retino_gamma_02','retino_movie','retino_rest'};
mask = load_untouch_nii('c:/shared/epireg/meanmask.nii.gz'); 
mean_nii = load_untouch_nii('c:/shared/epireg/mean.nii.gz'); 

maskinds = find(mask.img==1); 
for sb=1:length(subs)    
    corrs = load_untouch_nii(['e:/fmris/badger_',subs{sb},'/atlas_gamma_mcorrs.nii.gz']);
    corrs.img(isnan(corrs.img)) = 0; 
    [sv,si] = sort(corrs.img(:),'descend'); 
    for srcf=1:4
        
        cd(['E:/badger_eeg/',subs{sb}]);
        motion_si = load(['si_',eeg_set_names{srcf},'.mat']);
        motion_si = motion_si.si; 
        
        
        cd(['e:/fmris/badger_',subs{sb},'/',source_folders{srcf}]);        
        gamma_eeg = load_untouch_nii('gamma_atlas_hz.nii.gz');                
        res_gamma = reshape(gamma_eeg.img,[numel(gamma_eeg.img(:,:,:,1)),size(gamma_eeg.img,4)]); 
        mres_gamma = zscore(squeeze(mean(res_gamma(si(1:200),:),1))); 
        
        alpha_eeg = load_untouch_nii('alphabeta_atlas_hz.nii.gz');                
        res_alpha = reshape(alpha_eeg.img,[numel(alpha_eeg.img(:,:,:,1)),size(alpha_eeg.img,4)]); 
        mres_alpha = zscore(squeeze(mean(res_alpha(si(1:200),:),1))); 
        
        
        % motion correct gamma and alpha
        mc_mres_gamma = mres_gamma;
        mc_mres_alpha = mres_alpha; 
        
        mc_mres_gamma(motion_si(1:round(length(mres_gamma)*0.1))) = mean(mres_gamma);
        mc_mres_alpha(motion_si(1:round(length(mres_gamma)*0.1))) = mean(mres_alpha); 
                 
        sc_name = clean_fmri_names{srcf};        
        fmri = load_untouch_nii(['../atlas_fmri/',sc_name,'.nii.gz']);
        res_fmri = reshape(fmri.img,[numel(fmri.img(:,:,:,1)),size(fmri.img,4)]);
        fmri_vox = res_fmri(maskinds,:); 
        
        mres_bold = zscore(squeeze(mean(res_fmri(si(1:200),:),1))); 
        
        conv_alpha = conv(mres_alpha,spm_hrf(0.693),'full');
        conv_alpha = conv_alpha(1:length(mres_alpha)); 
        
        mc_conv_alpha = conv(mc_mres_alpha,spm_hrf(0.693),'full');
        mc_conv_alpha = mc_conv_alpha(1:length(mc_mres_alpha)); 
        
        conv_gamma = conv(mres_gamma,spm_hrf(0.693),'full');
        conv_gamma = conv_gamma(1:length(mres_gamma)); 
        
        mc_conv_gamma = conv(mc_mres_gamma,spm_hrf(0.693),'full');
        mc_conv_gamma = mc_conv_gamma(1:length(mc_mres_gamma)); 
        
        gamma_corrs(sb,srcf) = corr2(conv_gamma(20:end-20),mres_bold(20:length(conv_gamma)-20)); 
        alpha_corrs(sb,srcf) = corr2(conv_alpha(20:end-20),mres_bold(20:length(conv_alpha)-20)); 

        mc_gamma_corrs(sb,srcf) = corr2(mc_conv_gamma(20:end-20),mres_bold(20:length(mc_conv_gamma)-20)); 
        mc_alpha_corrs(sb,srcf) = corr2(mc_conv_alpha(20:end-20),mres_bold(20:length(mc_conv_alpha)-20)); 
    end

end

mgamma = mean(gamma_corrs(:,1:2),2); 
malpha = mean(alpha_corrs(:,1:2),2); 
mc_mgamma = mean(mc_gamma_corrs(:,1:2),2); 
mc_malpha = mean(mc_alpha_corrs(:,1:2),2); 

cgamma_corrs(:,1) = mgamma; cgamma_corrs(:,2:3) = gamma_corrs(:,3:4); 
calpha_corrs(:,1) = malpha; calpha_corrs(:,2:3) = alpha_corrs(:,3:4); 
mc_cgamma_corrs(:,1) = mc_mgamma; mc_cgamma_corrs(:,2:3) = mc_gamma_corrs(:,3:4); 
mc_calpha_corrs(:,1) = mc_malpha; mc_calpha_corrs(:,2:3) = mc_alpha_corrs(:,3:4); 

both_gamma(1,:,:) = cgamma_corrs; both_gamma(2,:,:) = mc_cgamma_corrs;
both_alpha(1,:,:) = calpha_corrs; both_alpha(2,:,:) = mc_calpha_corrs; 
stimnames = {'gamma','movie','rest'};

subplot(1,2,1);
barwitherr(squeeze(std(both_gamma,0,2))'./sqrt(8),squeeze(mean(both_gamma,2))'); hold on ;

subplot(1,2,2); 
barwitherr(squeeze(std(both_alpha,0,2))'./sqrt(8),squeeze(mean(both_alpha,2))'); hold on ;








