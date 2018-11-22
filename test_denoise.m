subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie'} ; 

for sub=1%:length(subs)
    cd(['c:/shared/badger_mri/',subs{sub},'/nii']) ; 
    ls
    hrf = load('gamma_hrf_1.txt') ; 
    
    bp = load_untouch_nii('bp_reg_topup_mc_retino_gamma_01.nii.gz') ; 
    reg = load_untouch_nii('reg_topup_mc_retino_gamma_01.nii.gz') ; 
    denoised = load_untouch_nii('allmelodic/denoised_reg_topup_mc_retino_gamma_01.nii.gz') ; 
    denoised_bp = load_untouch_nii('allmelodic/bp_denoised.nii.gz') ; 

    corrbp = voxcorr(squeeze(bp.img(:,:,:,10:end-10)),hrf(10:end-10)) ; 
    corrreg = voxcorr(squeeze(reg.img(:,:,:,10:end-10)),hrf(10:end-10)) ; 
    corrdenoised = voxcorr(squeeze(denoised.img(:,:,:,10:end-10)),hrf(10:end-10)) ; 
    corrdenoisedbp = voxcorr(squeeze(denoised_bp.img(:,:,:,10:end-10)),hrf(10:end-10)) ; 
    
    
    allts = load('reg_topup_mc_retino_gamma_01.txt') ; 
    
    bp_allts = eegfiltfft(allts',1/.693,.02,1) ; 
    c = corr(bp_allts(:,10:end-10)',hrf(10:end-10)) ;
    
    plot(mat2gray(bp_allts(64,10:end-10))) ; hold on ; plot(mat2gray(hrf(10:end-10)),'r') ; 
    
    
end
