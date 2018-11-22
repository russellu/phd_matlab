clear all ; close all ; 

subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 

for sb=1%:length(subs)
    
    cd(['E:\badger_eeg\',subs{sb}]);
    
    ls('*set'); 
    
    eeg = pop_loadset('retino_rest.set');
    trigs = {eeg.urevent.type};
    trtrigs = find(strcmpi('R128',trigs)); 
    lats = {eeg.urevent.latency}; 
    trlats = cell2mat(lats(trtrigs)); 

    freqs = 1:2:100; 
    
    for elec=1:64

        fcount = 1; 
        filtdat = zeros(50,size(eeg.data,2)); 
        for f=1:length(freqs)
            filtdat(fcount,:) = eegfiltfft(eeg.data(elec,:),eeg.srate,freqs(f)-3,freqs(f)+3);    
            fcount = fcount + 1; 
        end
        filtdat = filtdat(:,trlats(1):trlats(end)); 
        filtdat = imresize(filtdat,[50,450]); 
        allfilts(elec,:,:) = filtdat; 
    end
    
    
    logabsfilts = log(abs(allfilts)); 
    absfilts = (abs(allfilts)); 

    cd(['E:\fmris\badger_',subs{sb}]);
    fmri = load_untouch_nii('bp_reg_topup_mc_retino_rest.nii.gz'); 
    mask = load_untouch_nii('mask.nii.gz'); 
    maskinds = find(mask.img==1); 
    [mx,my,mz] = ind2sub(size(mask.img),maskinds); 
    
    globsig = mean(resimg,1);   
    clear allxcorrs; 
    for i=1:64
        for j=1:50
            allxcorrs(i,j,:) = xcorr(squeeze(absfilts(i,j,20:end-20)),globsig(20:end-20),20,'coeff'); % fmri should always go 2nd
        end
    end
    
    mabsfilts = squeeze(mean(absfilts,1)); 
    voxcorrs = zeros(length(maskinds),50,41); 
    for i=1:length(maskinds) ; disp(i); 
        for j=1:50
           ts = squeeze(fmri.img(mx(i),my(i),mz(i),:)); 
           voxcorrs(i,j,:) = xcorr(squeeze(mabsfilts(j,20:end-20)),ts(20:end-20)',20,'coeff');  
        end
    end
    
    mvoxcorrs = squeeze(mean(mean(voxcorrs(:,3:6,24:30),2),3)); 
    cmask = zeros(size(mask.img)); 
    for i=1:length(mvoxcorrs)
       cmask(mx(i),my(i),mz(i)) = mvoxcorrs(i);  
    end
    
    hrf = spm_hrf(1); 
    ideal = zeros(1,length(hrf)); ideal(1) = 1 ;
    xc = xcorr(hrf,ideal,40,'coeff'); 

        
end
    