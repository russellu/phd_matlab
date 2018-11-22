clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 
source_folders = {'den_retino_allstims_01','den_retino_allstims_02','den_retino_gamma_01','den_retino_gamma_02','den_retino_movie','den_retino_rest'};
clean_fmri_names = {'bp_clean_retino_allstims_01','bp_clean_retino_allstims_02','bp_clean_retino_gamma_01','bp_clean_retino_gamma_02','bp_clean_retino_movie','bp_clean_retino_rest'};
gsr_fmri_names = {'bp_gsr_clean_retino_allstims_01','bp_gsr_clean_retino_allstims_02','bp_gsr_clean_retino_gamma_01','bp_gsr_clean_retino_gamma_02','bp_gsr_clean_retino_movie','bp_gsr_clean_retino_rest'};

for sb=1:length(subs)
    
    mask = load_untouch_nii('c:/shared/epireg/meanmask.nii.gz'); 
    maskinds = find(mask.img==1); 
    ref50 = load_untouch_nii('e:/meanepis/ref50.nii.gz'); 
    
    
    for srcf=1:length(source_folders)
        cd(['e:/fmris/badger_',subs{sb},'/',source_folders{srcf}]);
        hzs = dir('atlas_hz*'); 
        
        sc_name = gsr_fmri_names{srcf};
        
        fmri = load_untouch_nii(['../atlas_fmri/',sc_name,'.nii.gz']);
        res_fmri = reshape(fmri.img,[numel(fmri.img(:,:,:,1)),size(fmri.img,4)]);
        fmri_vox = res_fmri(maskinds,:); 
        hz_corrs = zeros(size(fmri_vox,1),50); 
        for hz=1:length(hzs) ; disp([subs{sb},' ',num2str(hz)]); 
            eeg_src = load_untouch_nii(hzs(hz).name); 
            res_eeg = reshape(eeg_src.img,[numel(eeg_src.img(:,:,:,1)),size(eeg_src.img,4)]); 
            eeg_vox = res_eeg(maskinds,:); 
            conved_vox = zeros(size(eeg_vox)); 
            % hrf corrs
            hrf = spm_hrf(0.693); 
            
            for i=1:size(eeg_vox,1)
                conved_i = conv(eeg_vox(i,:),hrf,'full');
                conved_i = conved_i(1:size(conved_vox,2)); 
                conved_vox(i,:) = conved_i; 
                hz_corrs(i,hz) = corr2(conved_vox(i,20:size(eeg_vox,2)-20),fmri_vox(i,20:size(eeg_vox,2)-20)); 
            end
        end
        
        hz_3dimg = zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3),50); 
        [cx,cy,cz] = ind2sub(size(mask.img),maskinds); 
        for i=1:length(cx)
           hz_3dimg(cx(i),cy(i),cz(i),:) = hz_corrs(i,:);  
        end
        ref50.img = hz_3dimg; 
        save_untouch_nii(ref50,['hrfcorrs_',sc_name,'.nii.gz']); 
        
        
    end
    
end

%{
mask = load_untouch_nii('c:/shared/epireg/meanmask.nii.gz'); 
maskinds = find(mask.img==1); 
ref50 = load_untouch_nii('e:/meanepis/ref50.nii.gz'); 
corr_img = zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3),50,8,6); 
for sb=1:length(subs)
    
    for sf=1:length(source_folders)
        cd(['e:/fmris/badger_',subs{sb},'/',source_folders{sf}]); 
        hrf_name = dir('hrf*bp*gsr*');
        hrf_nii = load_untouch_nii(hrf_name.name); 
        corr_img(:,:,:,:,sb,sf) = hrf_nii.img; 
    end
    
end
corr_img(isnan(corr_img(:))) = 0 ;

cd E:\meanepis

ref50 = load_untouch_nii('ref50.nii.gz');
ref50.img = squeeze(mean(mean(corr_img,5),6)); 
save_untouch_nii(ref50,'mean_all_vox_gsr.nii.gz'); 


mcorrs = squeeze(mean(mean(corr_img(:,:,:,:,:,6),5),6)); 

for i=1:50 ; subplot(5,10,i) ; imagesc(imrotate(squeeze(mean(mcorrs(:,:,8:15,i),3)),270),[-.15,.15]) ; colormap jet; end
%}













