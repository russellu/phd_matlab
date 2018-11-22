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
                gamma_epochs(sb,i,j+16*(srcf-1),:) = mres_gamma(stimlats(i,j)-3:stimlats(i,j)+18); 
                alpha_epochs(sb,i,j+16*(srcf-1),:) = mres_alpha(stimlats(i,j)-3:stimlats(i,j)+18); 
                bold_epochs(sb,i,j+16*(srcf-1),:) = mres_bold(stimlats(i,j)-3:stimlats(i,j)+18); 
            end
        end           
    end

end

clear good_gamms good_alphs
for i=1:8
    for j=1:3
        good_gammas(i,j,:) = abs(zscore(mean(gamma_epochs(i,j,:,6:9),4))) <4; 
        good_alphas(i,j,:) = abs(zscore(mean(alpha_epochs(i,j,:,6:9),4))) <4; 
    end
end

clear gamma_corrs alpha_corrs
for i=1:8
    for j=1:3
           gamma_corrs(i,j) = corr2(squeeze(mean(bold_epochs(i,j,squeeze(good_gammas(i,j,:)),15:17),4)),squeeze(mean(gamma_epochs(i,j,squeeze(good_gammas(i,j,:)),6:12),4)));  
           alpha_corrs(i,j) = corr2(squeeze(mean(bold_epochs(i,j,squeeze(good_alphas(i,j,:)),15:17),4)),squeeze(mean(alpha_epochs(i,j,squeeze(good_alphas(i,j,:)),6:12),4)));  
    end
end


goodsubs = [1,2,3,5,6,7,8]; 
stimlabs = {'0%rnd','10%rnd','100%rnd'};

subplot(2,2,1);
barwitherr(std(gamma_corrs(goodsubs,:),0,1)/sqrt(8),mean(gamma_corrs(goodsubs,:)),'r'); title('gamma vs BOLD single trials'); ylabel('correlation (r)'); set(gca,'XTickLabel',stimlabs); 
subplot(2,2,2); 
barwitherr(std(alpha_corrs(goodsubs,:),0,1)/sqrt(8),mean(alpha_corrs(goodsubs,:))); title('alpha/beta vs BOLD single trials'); ylabel('correlation (r)');set(gca,'XTickLabel',stimlabs); 
subplot(2,2,3); 

res_alpha = reshape(alpha_corrs(goodsubs,:),[1,21]);
res_gamma = reshape(gamma_corrs(goodsubs,:),[1,21]); 
both_freqs = [res_gamma;res_alpha]; 
bar(1,mean(both_freqs(1,:),2),'r'); hold on ; bar(2,mean(both_freqs(2,:),2)); 
errorbar(mean(both_freqs,2),std(both_freqs,0,2)/sqrt(21),'k.'); 
title('mean all stim types'); ylabel('correlation(r)') ; set(gca,'XTick',[1,2],'XTickLabel',{'gamma','alpha/beta'});











%gamma_epochs = gamma_epochs - repmat(mean(gamma_epochs(:,:,:,1:3),4),[1,1,1,22]); 
%alpha_epochs = alpha_epochs - repmat(mean(alpha_epochs(:,:,:,1:3),4),[1,1,1,22]); 
%bold_epochs = bold_epochs - repmat(mean(bold_epochs(:,:,:,1:3),4),[1,1,1,22]); 










