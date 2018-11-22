clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 
source_folders = {'den_retino_gamma_01','den_retino_gamma_02','den_retino_movie','den_retino_rest'}; %'den_retino_allstims_01','den_retino_allstims_02',
clean_fmri_names = {'bp_clean_retino_gamma_01','bp_clean_retino_gamma_02','bp_clean_retino_movie','bp_clean_retino_rest'}; %'bp_clean_retino_allstims_01','bp_clean_retino_allstims_02',
gsr_fmri_names = {'bp_gsr_clean_retino_gamma_01','bp_gsr_clean_retino_gamma_02','bp_gsr_clean_retino_movie','bp_gsr_clean_retino_rest'};%'bp_gsr_clean_retino_allstims_01','bp_gsr_clean_retino_allstims_02',
eeg_set_names = {'retino_gamma_01','retino_gamma_02','retino_movie','retino_rest'};
mask = load_untouch_nii('c:/shared/epireg/meanmask.nii.gz'); 
maskinds = find(mask.img==1); 
for sb=1:length(subs)    
    
    corrs = load_untouch_nii(['e:/fmris/badger_',subs{sb},'/atlas_gamma_mcorrs.nii.gz']);
    corrs.img(isnan(corrs.img)) = 0; 
    [sv,si] = sort(corrs.img(:),'descend'); 

    for srcf=1:length(source_folders)
        cd(['e:/fmris/badger_',subs{sb},'/',source_folders{srcf}]);        
        gamma_eeg = load_untouch_nii('gamma_atlas_hz.nii.gz');        
        
        res_gamma = reshape(gamma_eeg.img,[numel(gamma_eeg.img(:,:,:,1)),size(gamma_eeg.img,4)]); 
        mres_gamma = zscore(squeeze(mean(res_gamma(si(1:200),:),1))); 
        
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
        
        if srcf <= 2                        
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
                    gamma_epochs(sb,i,j+16*(srcf-1),:) = mres_gamma(stimlats(i,j)-3:stimlats(i,j)+14); 
                    bold_epochs(sb,i,j+16*(srcf-1),:) = mres_bold(stimlats(i,j)-3:stimlats(i,j)+14); 
                end
            end        
        
        elseif srcf ==3
            movie_xcorrs(sb,:) = xcorr(mres_bold(20:end-20),mres_gamma(20:end-20),12,'coeff'); 
            
        elseif srcf ==4
            rest_xcorrs(sb,:) = xcorr(mres_bold(20:end-20),mres_gamma(20:end-20),12,'coeff'); 
        end
        %corrs(:,:,:,srcf) = voxcorr(fmri.img(:,:,:,20:length(trlats)-20),conved(20:length(trlats)-20)); 
    end
    %{
    gamma_epochs = gamma_epochs - repmat(mean(gamma_epochs(:,:,1:3),3),[1,1,18]); 
    bold_epochs = bold_epochs - repmat(mean(bold_epochs(:,:,1:8),3),[1,1,22]); 

    figure,
    colors = {'r','g','b'};
    for i=1:3        
       shadedErrorBar([],squeeze(mean(gamma_epochs(i,:,:),2)),squeeze(std(gamma_epochs(i,:,:),0,2))/sqrt(32),{colors{i}}) ; hold on; 
       shadedErrorBar([],squeeze(mean(bold_epochs(i,:,:),2)),squeeze(std(bold_epochs(i,:,:),0,2))/sqrt(32),{colors{i}}) ; hold on; 
    end
    
    %}
    %{
    mcorrs = squeeze(mean(corrs,4));     
    cd(['e:/fmris/badger_',subs{sb}]);        
    mask.img = mcorrs; 
    save_untouch_nii(mask,'atlas_gamma_mcorrs.nii.gz'); 
    %}
end

clear xcorrs; 
for i=1:8
    for j=1:3
        xcorrs(i,j,:) = xcorr(squeeze(mean(bold_epochs(i,j,:,:),3)),squeeze(mean(gamma_epochs(i,j,:,:),3)),12,'coeff'); 
    end
end

plot(squeeze(mean(xcorrs,1))') ; hline(0,'k') ; hold on ; 
hrf = spm_hrf(0.693); zhrf = zeros(1,size(xcorrs,3)); 
zhrf(13:end) = hrf(1:13); plot(zhrf*2,'k','LineWidth',2)
plot(mean(movie_xcorrs,1),'g','LineWidth',2);
plot(mean(rest_xcorrs,1),'b','LineWidth',2);





%{
gamma_epochs = gamma_epochs - repmat(mean(gamma_epochs(:,:,:,1:3),4),[1,1,1,18]); 
bold_epochs = bold_epochs - repmat(mean(bold_epochs(:,:,:,1:8),4),[1,1,1,22]); 
colors = {'r','g','b'};
subplot(1,2,1) ;
for i=1:3
    shadedErrorBar([],squeeze(mean(mean(gamma_epochs(:,i,:,:),1),3)),squeeze(std(mean(gamma_epochs(:,i,:,:),3),0,1))/sqrt(8),{colors{i}}); hold on ;     
    vline([4,11],'k'); 
end
title('Gamma (n=8)'); 
subplot(1,2,2); 
for i=1:3
    shadedErrorBar([],squeeze(mean(mean(bold_epochs(:,i,:,:),1),3)),squeeze(std(mean(bold_epochs(:,i,:,:),3),0,1))/sqrt(8),{colors{i}}); hold on ;    
    vline([4,11],'k'); 
end
title('BOLD (n=8)'); 
%}


