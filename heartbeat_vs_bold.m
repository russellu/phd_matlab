clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 
eeg_set_names = {'retino_gamma_01','retino_gamma_02','retino_movie','retino_rest'};
fmri_names = {'bp_clean_retino_gamma_01','bp_clean_retino_gamma_02','bp_clean_retino_movie','bp_clean_retino_rest'}; 
mask = load_untouch_nii('c:/shared/epireg/meanmask.nii.gz'); 
maskinds = find(mask.img==1); 

for sb=1:length(subs)    
    for set=1:length(eeg_set_names)
        cd(['E:/badger_eeg/',subs{sb}]); ls
        eeg = pop_loadset([eeg_set_names{set},'.set']);
        trigs = {eeg.urevent.type};
        trtrigs = find(strcmpi('R128',trigs)); 
        lats = {eeg.urevent.latency}; 
        trlats = cell2mat(lats(trtrigs)); 
        
        freqs = 1:2:100;
        clear resfilts
        for f=1:length(freqs)
           filts = eegfiltfft(eeg.data([1:31,33:64],:),eeg.srate,freqs(f)-1,freqs(f)+1);  
           filts = filts(:,trlats(1):trlats(end)); 
           resfilts(:,f,:) = imresize(abs(filts),[63,length(trlats)]); 
        end
        resfilts = squeeze(mean(resfilts,1)); 
        alpha_resfilts = squeeze(mean(resfilts(4:12,:),1)); 
        alpha_resfilts = eegfiltfft(alpha_resfilts,1/0.693,0.005,1); 
        
        bp_resfilts = eegfiltfft(resfilts,1/0.693,0.006,1); 
        clear hrf_bp_resfilts
        for i=1:size(bp_resfilts,1)
            conved = conv(bp_resfilts(i,:),spm_hrf(0.693),'full');
            hrf_bp_resfilts(i,:) = conved(1:size(bp_resfilts,2));             
        end
        
        
        cd(['E:/fmris/badger_',subs{sb},'/atlas_fmri']); ls 
        fmri = load_untouch_nii([fmri_names{set},'.nii.gz']);
        res_fmri = reshape(fmri.img,[numel(fmri.img(:,:,:,1)),size(fmri.img,4)]);
        masked_fmri = res_fmri(maskinds,:); 
        clear xcorrs
        for i=1:length(maskinds)
            xcorrs(i,:) = xcorr(alpha_resfilts(20:end-20),masked_fmri(i,20:size(resfilts,2)-20),20,'coeff');             
        end
        res_xcorrs = zeros(size(res_fmri,1),41);
        res_xcorrs(maskinds,:) = xcorrs; 
        res_xcorrs = reshape(res_xcorrs,[size(fmri.img,1),size(fmri.img,2),size(fmri.img,3),41]); 
        
        all_res_xcorrs(:,:,:,:,sb,set) = res_xcorrs; 
        
        clear hrfcorrs; 
        for j=1:50
            hrfcorrs(:,j) = corr(masked_fmri(:,20:size(hrf_bp_resfilts,2)-20)',hrf_bp_resfilts(j,20:end-20)'); 
        end
        res_hrfcorrs = zeros(size(res_fmri,1),50); 
        res_hrfcorrs(maskinds,:) = hrfcorrs; 
        res_hrfcorrs = reshape(res_hrfcorrs,[size(fmri.img,1),size(fmri.img,2),size(fmri.img,3),50]); 
              
        
        all_res_hrf_corrs(:,:,:,:,sb,set) = res_hrfcorrs; 
        
    end
    
end

all_res_hrf_corrs(isnan(all_res_hrf_corrs))= 0; 
mean_hrf = squeeze(mean(mean(all_res_hrf_corrs,5),6)); 
cd e:/meanepis ; ls 
ref50 = load_untouch_nii('ref50.nii.gz');
ref50.img = mean_hrf ; 
save_untouch_nii(ref50,'mean_bcg_corrs.nii.gz'); 


all_res_xcorrs(isnan(all_res_xcorrs)) = 0;
ref41=load_untouch_nii('ref41.nii.gz');
ref41.img = squeeze(mean(mean(all_res_xcorrs,5),6));
save_untouch_nii(ref41,'mean_bcg_xcorrs.nii.gz'); 






