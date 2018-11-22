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
      
        sc_name = clean_fmri_names{srcf};        
        fmri = load_untouch_nii(['../atlas_fmri/',sc_name,'.nii.gz']);
        res_fmri = reshape(fmri.img,[numel(fmri.img(:,:,:,1)),size(fmri.img,4)]);
        fmri_vox = res_fmri(maskinds,:); 
        
        mres_bold = zscore(squeeze(mean(res_fmri(si(1:200),:),1))); 
        
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
        
        bold_corrs(sb,srcf) = corr2(mres_bold(21:end-20),conv_zsummat(20:length(mres_bold)-21)); 

    end

end



