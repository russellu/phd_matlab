cd('C:\shared\eeg_badger\russ') ; ls ; clear all ; close all ; 
%filenames = {'Retino_angles_2.vhdr','Retino_angles_3.vhdr','Retino_angles_4.vhdr','Retino_angles_5.vhdr','Retino_angles_6.vhdr','Retino_polar_1.vhdr','rest.vhdr'} ; 

for fname = 1:length(filenames)
    EEG = pop_loadbv('.',filenames{fname}) ; 
    EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
    EEG = denoise_grad2(EEG) ; 
    %EEG = denoise_bcg(EEG) ; 
    eegs{fname} = EEG ; 
    if fname==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end   
end
mergefilt = merged ; mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,mergefilt.srate,0,0.5) ; 
mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,mergefilt.srate,59.5,60.5) ; 
ica = pop_runica(mergefilt,'runica') ;

% get the FMRI files:
cd c:/shared/badger/russ/mel
fmrinames = {'../bp_reg_mc_Test_Russell_2015_10_15_WIP_EEG-fMRI_MB2_3mm_SENSE_12_1.nii.gz'} ; 
%fmrinames = {'reg_mc_fmri_1.nii.gz','reg_mc_fmri_2.nii.gz','reg_mc_fmri_3.nii.gz','reg_mc_fmri_4.nii.gz','reg_mc_fmri_5.nii.gz','reg_mc_fmri_6.nii.gz','reg_mc_fmri_7.nii.gz',} ; 
for fname=1
    cfmri = load_untouch_nii(fmrinames{1}) ; fimg = cfmri.img ; 
    ceeg = eegs{fname} ; 
    cevents = {ceeg.urevent.type} ; 
    clats = {ceeg.urevent.latency} ; 
    r128s = find(strcmp('R128',cevents)) ; r128lats = round(cell2mat(clats(r128s))) ; 
    ceeg.icaact = icaact(ceeg.data,ica.icaweights*ica.icasphere,[]) ; 
    ceeg.icachansind = ica.icachansind ; ceeg.icawinv = ica.icawinv ; ceeg.icaweights = ica.icaweights ; ceeg.icasphere = ica.icasphere ; 
    
    difflats = diff(r128lats) ; trlen = round(difflats(1)) ; 
    icadata = ceeg.icaact(:,r128lats(1):r128lats(length(r128lats))+trlen) ;

    freqs = 1:80 ; 
    for f=1:length(freqs)-1
        data(f,:,:) = eegfiltfft(icadata,ica.srate,freqs(f),freqs(f+1)) ; 
    end
    
    sqrdata = sqrt(data.^2) ; 
    % get the power in increments of TR
    tr = cfmri.hdr.dime.pixdim(5) ; % TR in secs
    winlength = round(tr*ceeg.srate) ; 
    powdata = zeros(size(data,1),size(data,2),410) ; 
    for i=1:size(data,1)
        for j=1:size(data,2) ; kcount = 1 ; 
            for k=1:winlength:(size(data,3))-winlength
                powdata(i,j,kcount) = squeeze(mean(sqrdata(i,j,k:k+winlength),3)) ; 
                kcount = kcount + 1 ; 
            end
            powdata(i,j,kcount) = squeeze(mean(powdata(i,j,:),3)) ; 
        end
    end
       
    melmix = load('melodic_mix') ;    
    restest = imresize(squeeze(powdata(5,20,:)),[size(icadata,2),1]) ;  
    
    conversp = zeros(size(powdata)) ; hrf = spm_hrf((trlen/ceeg.srate)) ; 
    for i=1:size(powdata,1) 
        for j=1:size(powdata,2)
            convij = conv(squeeze(powdata(i,j,:)),hrf) ; 
            conversp(i,j,:) = convij(1:410) ; 
        end
    end
        
    mix1 = melmix(1:410,:) ; trimg = fimg(:,:,:,10:size(mix1,1)-20) ; mix1 = mix1(10:size(mix1,1)-20,:) ; 
    conversp = conversp(:,:,10:size(conversp,3)-20) ;
    for i=1:size(mix1,2)
       mix1(:,i) = detrend(mix1(:,i),'linear') ;  
    end
    % correlate each EEG component with each FMRI component at each frequency
    clear corrs ;
    for i=1:64 ;
       corrs(i,:,:) = corr(squeeze(conversp(:,i,:))',mix1) ;          
    end
    fcomp = 2 ; ecomp = 22 ; freqind = 14 ; 
    etrace = squeeze(conversp(ecomp,:,freqind)) ;    
    figure, plot(mat2gray(squeeze(conversp(ecomp,freqind,:))),'r','LineWidth',3) ; hold on ; plot(smooth(mat2gray(mix1(:,fcomp))),'b','LineWidth',3) ; 
    title(num2str(corr2(squeeze(conversp(ecomp,:,freqind)),mix1(:,fcomp)'))) ; 

    restrimg = reshape(trimg,[size(trimg,1)*size(trimg,2)*size(trimg,3),size(trimg,4)]) ; 
    detr = detrend(restrimg','linear')' ; detres = reshape(detr,size(trimg)) ; 
    %{
    for i=1:64 ; 
        etrace = squeeze(mean(conversp(5,i,:),1)) ; disp(i) ; 
        corrbrain = voxcorr(detres,etrace) ; meanf = load_untouch_nii('mean.nii.gz') ; 
        meanf.img = corrbrain ; save_untouch_nii(meanf,['corrbrain_,',num2str(i),'.nii.gz']) ; 
    end
    %}
    meanf = load_untouch_nii('mean.nii.gz') ; 
    for i=1:size(conversp,1) ; 
        etrace = squeeze(mean(conversp(i,23,:),1)) ; disp(i) ; 
        corrbrain = voxcorr(detres,etrace) ; 
        meanf.img = corrbrain ; save_untouch_nii(meanf,['corrbrain_,',num2str(i),'.nii.gz']) ; 
    end
end




