clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','sukhman','valerie'}; 
scans = {'bp_reg_topup_mc_retino_allstims_01','bp_reg_topup_mc_retino_allstims_02','bp_reg_topup_mc_retino_gamma_01','bp_reg_topup_mc_retino_gamma_02','bp_reg_topup_mc_retino_movie','bp_reg_topup_mc_retino_rest'}; 
eegscans = {'retino*allstim*01*set','retino*allstim*02*set','retino*gamma*01*set','retino*gamma*02*set','retino*movie*set','retino*rest*set'}; 

for sb=1:length(subs)
    ls 
    %figure,
    for sc=1:length(scans)
            cd(['E:\rawbadger\badger_mri\',subs{sb},'\nii\warped']);   

        nii = load_untouch_nii([scans{sc},'.nii.gz']);  

        freqs = abs(fft(nii.img,[],4)); 
        freqs = freqs(:,:,:,1:end/2-1); 
        mfreqs = squeeze(mean(freqs(:,:,:,120:220),4)); 
        allfreqs(sb,sc,:,:,:) = mfreqs; 


        [sv,si] = sort(mfreqs(:),'descend'); 
        res_nii = reshape(nii.img,[numel(nii.img(:,:,:,1)),size(nii.img,4)]); 
        % subplot(2,3,sc);
        % plot(squeeze(mean(res_nii(si(1:200),:)))); 

        bold_hr = squeeze(mean(res_nii(si(1:200),:))); 
        
        cd(['E:\badger_eeg\',subs{sb}]);  
        eegscan = dir(eegscans{sc}); 
        disp(eegscan.name); 

        eeg = pop_loadset(eegscan.name); 
        raweeg = eeg; 
        trigs = {eeg.urevent.type};
        trtrigs = find(strcmpi('R128',trigs)); 
        lats = {eeg.urevent.latency}; 
        trlats = cell2mat(lats(trtrigs)); 

        freqs = 1:2:100; 
        clear res_filts bcg_res_filts    
        for f=1:length(freqs)
           filt_f = eegfiltfft(eeg.data,eeg.srate,freqs(f)-1,freqs(f)+1);  
           filt_f = filt_f(:,trlats(1):trlats(end)); 
           res_filts(f,:,:) = imresize(abs(filt_f),[64,length(trlats)]);  
        end

        for i=1:50
           res_filts(i,:,:) = eegfiltfft(squeeze(res_filts(i,:,:)),1/0.693,0.005,2);  
        end

        
        meanmask = load_untouch_nii('c:/shared/epireg/meanmask.nii.gz'); 
        maskind = find(meanmask.img==1); 
        globsig = mean(res_nii(maskind,:)); 
        
        globsig = globsig(1:length(trlats)); 
        bold_hr = bold_hr(1:length(trlats)); 
        
        for i=1:50
            for j=1:64
                hr_xcorrs(sb,sc,i,j,:) = xcorr(bold_hr(20:end-20),squeeze(res_filts(i,j,20:end-20)),20,'coeff'); 
                global_xcorrs(sb,sc,i,j,:) = xcorr(globsig(20:end-20),squeeze(res_filts(i,j,20:end-20)),20,'coeff'); 
            end
        end        

        
        
        
        
        
        
    end
    
    
    
    
    
    
end







