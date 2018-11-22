clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 
scans = {'bp_reg_topup_mc_retino_allstims_01','bp_reg_topup_mc_retino_allstims_02','bp_reg_topup_mc_retino_gamma_01','bp_reg_topup_mc_retino_gamma_02','bp_reg_topup_mc_retino_movie','bp_reg_topup_mc_retino_rest'}; 
eegscans = {'retino*allstim*01*set','retino*allstim*02*set','retino*gamma*01*set','retino*gamma*02*set','retino*movie*set','retino*rest*set'}; 

for sb=1%:length(subs)
    
    for sc=1%:length(scans)
        cd(['E:\badger_eeg\',subs{sb}]);  
        eegscan = dir(eegscans{sc}); 
        disp(eegscan.name); 
        grad = pop_loadset(eegscan.name); 
        
        bcg_chan = grad.data(32,:); 
        grad.data(32,:) = rand(1,size(grad.data,2))*0.001 ; 
        
        filt_bcg = smooth(abs(eegfiltfft(bcg_chan,grad.srate,3,30)),100); 
        
        % epoch the BCG
        bdat = grad.data ; bdat = eegfiltfft(bdat,250,0.5,125) ;
        nbdat = eegfiltfft(grad.data,250,3,125) ; 
        subdat = zeros(size(bdat)) ; 
        [weights,sphere] = runica(bdat,'maxsteps',128) ; 
        acts = weights*sphere*bdat ; 

        filt_acts = imfilter(abs(eegfiltfft(acts,grad.srate,3,30)),fspecial('gaussian',[1,100],25)); 
        
        for i=1:10
            xcorrs(i,:) = xcorr(filt_acts(i,end/8:end-end/8),filt_bcg(end/8:end-end/8),200,'coeff'); 
        end
        max_xcorrs = max(xcorrs,[],2); 
        max_ind = find(max_xcorrs==max(max_xcorrs)); 
        
        
    end
    
    
end