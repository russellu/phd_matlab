cd('C:\shared\badger\russ\russ') ; ls ; clear all ; close all ; 
%filenames = {'Retino_angles_2.vhdr','Retino_angles_3.vhdr','Retino_angles_4.vhdr','Retino_angles_5.vhdr','Retino_angles_6.vhdr','Retino_polar_1.vhdr','rest.vhdr'} ; 
filenames=dir('*vhdr') ; 
for fname = 1:length(filenames)
    EEG = pop_loadbv('.',filenames(fname).name) ; 
    EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
    EEG = denoise_grad2(EEG) ; 
    EEG = denoise_bcg(EEG) ; 
    eegs{fname} = EEG ; 
    if fname==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end   
end
mergefilt = merged ; mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,mergefilt.srate,0,0.5) ; 
mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,mergefilt.srate,59.5,60.5) ; 
ica = pop_runica(mergefilt,'runica') ;

% get the FMRI files:
cd c:/shared/badger/russ
fmrinames = {'bp_mc_Test_Russell_2015_10_15_WIP_EEG-fMRI_MB2_3mm_SENSE_13_1.nii.gz'} ; 
%fmrinames = {'reg_mc_fmri_1.nii.gz','reg_mc_fmri_2.nii.gz','reg_mc_fmri_3.nii.gz','reg_mc_fmri_4.nii.gz','reg_mc_fmri_5.nii.gz','reg_mc_fmri_6.nii.gz','reg_mc_fmri_7.nii.gz',} ; 
for fname=2
    cfmri = load_untouch_nii(fmrinames{1}) ; fimg = cfmri.img ; 
    ceeg = eegs{fname} ; 
    cevents = {ceeg.urevent.type} ; 
    clats = {ceeg.urevent.latency} ; 
    r128s = find(strcmp('R128',cevents)) ; r128lats = round(cell2mat(clats(r128s))) ; 
    ceeg.icaact = icaact(ceeg.data,ica.icaweights*ica.icasphere,[]) ; 
    ceeg.icachansind = ica.icachansind ; ceeg.icawinv = ica.icawinv ; ceeg.icaweights = ica.icaweights ; ceeg.icasphere = ica.icasphere ; 
    
    difflats = diff(r128lats) ; trlen = round(difflats(1)) ; 
    icadata = ceeg.icaact(:,r128lats(1):r128lats(length(r128lats))+trlen) ;

    freqs = 1:2:80 ; clear data
    for f=1:length(freqs)-1
        data(f,:,:) = eegfiltfft(icadata,ica.srate,freqs(f),freqs(f+1)) ; 
    end
    
    sqrdata = (data.^2) ; 
    % get the power in increments of TR
    tr = cfmri.hdr.dime.pixdim(5) ; % TR in secs
    winlength = round(tr*ceeg.srate) ; 
    powdata = zeros(size(data,1),size(data,2),size(fimg,3)) ; 
    for i=1:size(data,1)
        for j=1:size(data,2) ; kcount = 1 ; 
            for k=1:winlength:(size(data,3))-winlength
                powdata(i,j,kcount) = squeeze(mean(sqrdata(i,j,k:k+winlength),3)) ; 
                kcount = kcount + 1 ; 
            end
            powdata(i,j,kcount) = squeeze(mean(powdata(i,j,:),3)) ; 
        end
    end
           
    conversp = zeros(size(powdata)) ; hrf = spm_hrf((trlen/ceeg.srate)) ; 
    for i=1:size(powdata,1) 
        for j=1:size(powdata,2)
            convij = conv(squeeze(powdata(i,j,:)),hrf) ; 
            conversp(i,j,:) = convij(1:size(powdata,3)) ; 
        end
    end
    
    conversp = conversp(:,:,10:size(powdata,3)-10) ;
    trimg = fimg(:,:,:,10:size(powdata,3)-10) ; %mix1 = mix1(10:size(powdata,3)-20,:) ;
    restrimg = reshape(trimg,[size(trimg,1)*size(trimg,2)*size(trimg,3),size(trimg,4)]) ; 
    detr = detrend(restrimg','linear')' ; detres = reshape(detr,size(trimg)) ; 
    
    meanf = load_untouch_nii('f1_mc_Test_Russell_2015_10_15_WIP_EEG-fMRI_MB2_3mm_SENSE_12_1.nii.gz') ; 
    for i=1:size(conversp,2) ; 
        etrace = squeeze(mean(conversp(1:3,i,:),1)) ; disp(i) ; 
        corrbrain = voxcorr(detres,etrace) ; 
        meanf.img = corrbrain ; save_untouch_nii(meanf,['corrbrain_',num2str(i),'.nii.gz']) ; 
    end
end




