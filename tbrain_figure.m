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
    clear gamma_epochs
    for srcf=1:2
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
        
        cd(['e:/eegs/badger_',subs{sb}]);
        eegname = dir([eeg_set_names{srcf},'*set']); 
        eeg = pop_loadset(eegname.name); 
        
        trigs = {eeg.urevent.type};
        trtrigs = find(strcmpi('R128',trigs)); 
        lats = {eeg.urevent.latency}; 
        trlats = cell2mat(lats(trtrigs)); 
        
        ts = zeros(1,size(eeg.data,2)); 
        
        stims = {'S  1','S  2','S  3'}; 
        for i=1:length(stims)
            stiminds_i = find(strcmpi(trigs,stims{i})); 
            lats_i = cell2mat(lats(stiminds_i)); 
            for j=1:length(lats_i)
               ts(lats_i(j):lats_i(j)+eeg.srate*5) = 1;  
               stimlats(i,j) = lats_i(j); 
            end
        end
        ts = ts(trlats(1):trlats(end)); 
        ts = imresize(ts,[1,length(trlats)]); 
        conved = conv(ts,spm_hrf(0.693),'full'); 
        conved = conved(1:length(ts)); 

        stimlats = stimlats - trlats(1); 
        stimlats = stimlats/eeg.srate; 
        stimlats = round(stimlats/0.693); 
        for i=1:size(stimlats,1)         
            for j=1:size(stimlats,2)
                gamma_epochs(i,j+16*(srcf-1),:,:) = res_gamma(:,stimlats(i,j)-3:stimlats(i,j)+18); 
                alpha_epochs(i,j+16*(srcf-1),:,:) = res_alpha(:,stimlats(i,j)-3:stimlats(i,j)+18); 
                bold_epochs(i,j+16*(srcf-1),:,:) = res_fmri(:,stimlats(i,j)-3:stimlats(i,j)+18); 
            end
        end        
        
       
    end
    
    
    gamma_epochs(isnan(gamma_epochs)) = 0 ;
    for i=1:size(gamma_epochs,1)
            [h,p,ci,stats] = ttest(squeeze(mean(gamma_epochs(i,:,:,5:12),4)),squeeze(mean(gamma_epochs(i,:,:,1:3),4))); 
            gammatbrains(i,:,:,:) = reshape(stats.tstat,size(corrs.img));
            [h,p,ci,stats] = ttest(squeeze(mean(alpha_epochs(i,:,:,5:12),4)),squeeze(mean(alpha_epochs(i,:,:,1:3),4))); 
            alphatbrains(i,:,:,:) = reshape(stats.tstat,size(corrs.img));            
            [h,p,ci,stats] = ttest(squeeze(mean(bold_epochs(i,:,:,15:18),4)),squeeze(mean(bold_epochs(i,:,:,5:10),4))); 
            boldtbrains(i,:,:,:) = reshape(stats.tstat,size(corrs.img));            
    end
    
    all_gammatbrains(sb,:,:,:,:) = gammatbrains; 
    all_alphatbrains(sb,:,:,:,:) = alphatbrains; 
    all_boldtbrains(sb,:,:,:,:) = boldtbrains; 

end


cd e:/tbrains ; 
figure,
all_gammatbrains(isnan(all_gammatbrains)) = 0; 
mtbrains = squeeze(mean(all_gammatbrains,1)); 
mtbrains(:,1,1,:) = 10; mtbrains(:,1,2,:) = -10; 
mtbrains(:,:,1,1) = 10; mtbrains(:,:,1,1) = -10; 
for i=1:3
    subplot(1,3,i)
    plotoverlayIntensity2D(squeeze(mean_nii.img(30,:,:)),squeeze(mat2gray(abs(mtbrains(i,30,:,:)))),squeeze(mtbrains(i,30,:,:)),90);
end
for i=1:3 ; mean_nii.img = permute(mtbrains,[2,3,4,1]);  mean_nii.img = mean_nii.img(:,:,:,i) ; save_untouch_nii(mean_nii,['gamma_t_',num2str(i),'.nii.gz']); end

figure,
all_alphatbrains(isnan(all_alphatbrains)) = 0; 
mtbrains = squeeze(mean(all_alphatbrains,1)); 
mtbrains(:,1,1,:) = 5; mtbrains(:,1,2,:) = -5; 
mtbrains(:,:,1,1) = 5; mtbrains(:,:,1,1) = -5; 
for i=1:3
    subplot(1,3,i)
    plotoverlayIntensity2D(squeeze(mean_nii.img(30,:,:)),squeeze(mat2gray(abs(mtbrains(i,30,:,:)))),squeeze(mtbrains(i,30,:,:)),90);
end
for i=1:3 ; mean_nii.img = permute(mtbrains,[2,3,4,1]);  mean_nii.img = mean_nii.img(:,:,:,i) ; save_untouch_nii(mean_nii,['alphabeta_t_',num2str(i),'.nii.gz']); end

figure,
all_boldtbrains(isnan(all_boldtbrains)) = 0; 
mtbrains = squeeze(mean(all_boldtbrains,1)); 
mtbrains(:,1,1,:) = 3; mtbrains(:,1,2,:) = -3; 
mtbrains(:,:,1,1) = 3; mtbrains(:,:,1,1) = -3; 
for i=1:3
    subplot(1,3,i)
    plotoverlayIntensity2D(squeeze(mean_nii.img(30,:,:)),squeeze(mat2gray(abs(mtbrains(i,30,:,:)))),squeeze(mtbrains(i,30,:,:)),90);
end
for i=1:3 ; mean_nii.img = permute(mtbrains,[2,3,4,1]); mean_nii.img = mean_nii.img(:,:,:,i) ; save_untouch_nii(mean_nii,['bold_t_',num2str(i),'.nii.gz']); end





