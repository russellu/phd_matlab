clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 
scans = {'bp_reg_topup_mc_retino_allstims_01','bp_reg_topup_mc_retino_allstims_02','bp_reg_topup_mc_retino_gamma_01','bp_reg_topup_mc_retino_gamma_02','bp_reg_topup_mc_retino_movie','bp_reg_topup_mc_retino_rest'}; 
shortnames = {'allstims_01','allstims_02','gamma_01','gamma_02','movie','rest'};
eegscans = {'retino*allstim*01*set','retino*allstim*02*set','retino*gamma*01*set','retino*gamma*02*set','retino*movie*set','retino*rest*set'}; 
comps = {[28,33,37],[18,46],[10,25,27],[11,32],[13,25],[28,39],[44,48,49],[35,46]};

for sb=1:length(subs)
    
    for scan=1:length(eegscans)
        cd(['E:\badger_eeg\',subs{sb}]);
        highweights = load('highweights.mat');
        highweights = highweights.highweights; 
        
        eeg_name = dir(eegscans{scan}); 
        
        gamma1 = pop_loadset(eeg_name.name); 

        winv = pinv(highweights{1}*highweights{2}); 
        acts = highweights{1}*highweights{2}*gamma1.data; 

        bads = zeros(1,64); bads(comps{sb}) = 1 ; bads = find(bads==0); 
        acts(bads,:) = 0;
        invacts = winv*acts; 

        trigs = {gamma1.urevent.type};
        trtrigs = find(strcmpi('R128',trigs)); 
        lats = {gamma1.urevent.latency}; 
        trlats = cell2mat(lats(trtrigs)); 

        clear all_filts; 
        freqs = 1:2:100;
        for f=1:length(freqs)
            filts_f = eegfiltfft(invacts,gamma1.srate,freqs(f)-1,freqs(f)+1); 
            res_f = imresize((abs(filts_f(:,round(trlats(1)):round(trlats(end))))),[64,length(trlats)]); 
            all_filts(f,:,:) = res_f; 
        end

        postelecs = [15,51,7,37,19,38,8,52,16,49,45,31,46,60,9,20,10];

        mpost = (log(squeeze(mean(all_filts(:,postelecs,:),2)))); 
        zmpost = zscore(mpost,[],2); 
        for i=1:size(zmpost,1)
           conved = conv(zmpost(i,:),spm_hrf(0.693),'full');
           zmpost(i,:) = conved(1:size(zmpost,2));
        end

        cd(['E:\rawbadger\badger_mri\',subs{sb},'\nii\warped']);   
        
        nii = load_untouch_nii([scans{scan},'.nii.gz']) ;
        nii.img = nii.img(:,:,:,1:length(trlats)); 
        clear voxcorrs
        for i=1:size(zmpost,1)
            voxcorrs(:,:,:,i) = voxcorr(nii.img(:,:,:,20:end-20),zmpost(i,20:end-20));  
        end
        ref50 = load_untouch_nii('E:\meanepis\ref50.nii.gz');
        ref50.img = voxcorrs;
        save_untouch_nii(ref50,['eegcorrs_',shortnames{scan},'.nii.gz']);
    end

end
