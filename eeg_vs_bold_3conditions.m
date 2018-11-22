subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie'} ; 

for sub=1%1:length(subs)

    %%%% BOLD FMRI processing:
    cd(['c:/shared/badger_mri/',subs{sub},'/nii']) ; 
    ls
    bold = load_untouch_nii('bp_reg_topup_mc_retino_gamma_01.nii.gz') ; 
    img = bold.img ; 
    
    %%%% EEG processing:
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    ls
    
    eegs = dir('1hz_preproc_*gamma*1*set') ;  
    EEG = pop_loadset(eegs(1).name) ; 
    %%% ground truth for each subject: epoch and check the gamma ERSP
    trigs = {'S  1','S  2','S  3'} ; 
    clear ersp
    for trig=1:3
        ep = pop_epoch(EEG,{trigs{trig}},[-2,7]) ; 
        for c=1:64 
            [ersp(trig,c,:,:),itc,powbase,times,freqs,~,~] = newtimef(ep.icaact(c,:),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                'plotersp','off','plotitc','off','baseline',0,'timesout',200,'nfreqs',120,'freqs',[1,120],'winsize',64) ; 
        end
    end
    for i=1:3 ; figure ;for j=1:64 ; subplot(5,13,j) ; imagesc(squeeze(ersp(i,j,:,:)),[-6,6]) ; title(j) ; end ; end
    %%% lots of noise in the data. need to remove noisy epochs in analyzer,
    %%% and then re-run ICA. especially on valerie
    trigs = {EEG.urevent.type} ; lats = cell2mat({EEG.urevent.latency}) ; 
    r128s = find(strcmp(trigs,'R128')) ;
    rlats = lats(r128s) ; 
    
    etrim = pop_select(EEG,'time',[rlats(1)./EEG.srate,rlats(length(rlats))./EEG.srate + EEG.srate*0.693]) ; 
    
    % get the continuous power in all frequencies
    for c=1:64 
        [epow(c,:,:),itc,powbase,times,freqs,~,~] = newtimef(etrim.icaact(c,:),etrim.pnts,[etrim.xmin,etrim.xmax],etrim.srate,0,...
            'plotersp','off','plotitc','off','baseline',NaN,'timesout',735,'nfreqs',120,'freqs',[1,120],'winsize',etrim.srate) ; 
    end
    imagesc(squeeze(mean(epow([19,22,31],:,:),1))); 
    
    gamma = zscore(squeeze(mean(epow(:,40:70,:),2))) ; 
    hrf = spm_hrf(.693) ; 
    for i=1:size(gamma,1)
        c = conv(squeeze(gamma(i,:)),hrf,'full') ; 
        conved(i,:) = detrend(c(1:735)) ;       
    end
    
    for i=1:64
        corrbrains(:,:,:,i) = voxcorr(img(:,:,:,10:end-10),conved(i,10:end-10)) ;      
    end
    save_nii(make_nii(corrbrains),'corrbrains.nii.gz') ; 
end







