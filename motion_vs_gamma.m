clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 
source_folders = {'den_retino_gamma_01','den_retino_gamma_02','den_retino_movie','den_retino_rest'}; %'den_retino_allstims_01','den_retino_allstims_02',
clean_fmri_names = {'bp_clean_retino_gamma_01','bp_clean_retino_gamma_02','bp_clean_retino_movie','bp_clean_retino_rest'}; %'bp_clean_retino_allstims_01','bp_clean_retino_allstims_02',
gsr_fmri_names = {'bp_gsr_clean_retino_gamma_01','bp_gsr_clean_retino_gamma_02','bp_gsr_clean_retino_movie','bp_gsr_clean_retino_rest'};%'bp_gsr_clean_retino_allstims_01','bp_gsr_clean_retino_allstims_02',
eeg_set_names = {'retino_gamma_01','retino_gamma_02','retino_movie','retino_rest'};
scans = {'mc_retino_gamma_01','mc_retino_gamma_02','mc_retino_movie','mc_retino_rest'}; 
mask = load_untouch_nii('c:/shared/epireg/meanmask.nii.gz'); 
mean_nii = load_untouch_nii('c:/shared/epireg/mean.nii.gz'); 

maskinds = find(mask.img==1); 
for sb=1:length(subs)    
    corrs = load_untouch_nii(['e:/fmris/badger_',subs{sb},'/atlas_gamma_mcorrs.nii.gz']);
    corrs.img(isnan(corrs.img)) = 0; 
    [sv,si] = sort(corrs.img(:),'descend'); 
    for srcf=1:4
        cd(['e:/fmris/badger_',subs{sb},'/',source_folders{srcf}]);        
        gamma_eeg = load_untouch_nii('gamma_atlas_hz.nii.gz');                
        res_gamma = reshape(gamma_eeg.img,[numel(gamma_eeg.img(:,:,:,1)),size(gamma_eeg.img,4)]); 
        mres_gamma = zscore(squeeze(mean(res_gamma(si(1:200),:),1))); 
        
        alpha_eeg = load_untouch_nii('alphabeta_atlas_hz.nii.gz');                
        res_alpha = reshape(alpha_eeg.img,[numel(alpha_eeg.img(:,:,:,1)),size(alpha_eeg.img,4)]); 
        mres_alpha = zscore(squeeze(mean(res_alpha(si(1:200),:),1))); 
        
        sc_name = clean_fmri_names{srcf};        
        fmri = load_untouch_nii(['../atlas_fmri/',sc_name,'.nii.gz']);
        res_fmri = reshape(fmri.img,[numel(fmri.img(:,:,:,1)),size(fmri.img,4)]);
        fmri_vox = res_fmri(maskinds,:); 
        
        mres_bold = zscore(squeeze(mean(res_fmri(si(1:200),:),1))); 
        
        conv_alpha = conv(mres_alpha,spm_hrf(0.693),'full');
        conv_alpha = conv_alpha(1:length(mres_alpha)); 
        
        conv_gamma = conv(mres_gamma,spm_hrf(0.693),'full');
        conv_gamma = conv_gamma(1:length(mres_gamma)); 
        
        %gamma_corrs(sb,srcf) = corr2(conv_gamma(20:end-20),mres_bold(20:length(conv_gamma)-20)); 
        %alpha_corrs(sb,srcf) = corr2(conv_alpha(20:end-20),mres_bold(20:length(conv_alpha)-20)); 
        
        
        cd(['E:\rawbadger\badger_mri\',subs{sb},'\nii\',scans{srcf},'.nii.mat']);
        mats = dir('MAT*');
        ls 
        
        for mat=1:length(mats)
            m1 = load(mats(mat).name);
            allmats(:,:,mat) = m1;    
        end

        clear newmat; 
        newmat(1,:) = squeeze(allmats(1,2,:)); 
        newmat(2,:) = squeeze(allmats(1,3,:)); 
        newmat(3,:) = squeeze(allmats(2,1,:)); 
        newmat(4,:) = squeeze(allmats(2,3,:)); 
        newmat(5,:) = squeeze(allmats(3,1,:)); 
        newmat(6,:) = squeeze(allmats(3,2,:)); 

        zmat = zscore(newmat,[],2); 

        summat = sum(abs(diff(zmat,1,2))); 
        zsummat = zscore(summat); 
        
        conv_zsummat = conv(zsummat,spm_hrf(0.693),'full'); 
        conv_zsummat = conv_zsummat(1:length(zsummat)); 
        
        gamma_corrs(sb,srcf) = corr2(conv_gamma(21:end-20),conv_zsummat(20:length(conv_gamma)-21)); 
        alpha_corrs(sb,srcf) = corr2(conv_alpha(21:end-20),conv_zsummat(20:length(conv_alpha)-21)); 

    end

end

mgamma = mean(gamma_corrs(:,1:2),2); 
malpha = mean(alpha_corrs(:,1:2),2); 

cgamma_corrs(:,1) = mgamma; cgamma_corrs(:,2:3) = gamma_corrs(:,3:4); 
calpha_corrs(:,1) = malpha; calpha_corrs(:,2:3) = alpha_corrs(:,3:4); 

stimnames = {'gamma','movie','rest'};
subplot(2,2,1);
barwitherr(squeeze(std(cgamma_corrs,0,1))/sqrt(8),mean(cgamma_corrs,1),'r'); set(gca,'XTickLabel',stimnames); ylabel('correlation (r)'); title('gamma vs motion continuous scan');
subplot(2,2,2); 
barwitherr(squeeze(std(calpha_corrs,0,1))/sqrt(8),mean(calpha_corrs,1)); set(gca,'XTickLabel',stimnames); ylabel('correlation (r)'); title('alpha vs motion continuous scan');

subplot(2,2,3); 
res_alpha = (calpha_corrs(:)); 
res_gamma = (cgamma_corrs(:)); 
res_both = [res_gamma,res_alpha]; 
bar(1,mean(res_gamma),'r') ; hold on ; bar(2,mean(res_alpha)); 
errorbar(mean(res_both,1),std(res_both,0,1)/sqrt(24),'.k'); set(gca,'XTick',[1,2],'XTickLabel',{'gamma','alpha/beta'}); title('mean all scans'); ylabel('correlation (r)'); 








