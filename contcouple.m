clear all ; close all ; 
fmrivols = {'blur_reg_bp_mc_retino_allstims_01.nii.gz','blur_reg_bp_mc_retino_allstims_02.nii.gz'} ; 
stimtrigs = {'S  1','S  2','S  3'} ; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie'} ;
subcomps = {[28,24,20,32],[35,39],[22,23,52],[25,44],[27,26,29],[19,26,21],[44,21]} ;
allcomps = {[17,21,39]} ; 
hrf = spm_hrf(0.693) ; 
for subby=1 %:7 ;
    clear minvals mininds corrvols diffvols taskvols basevols allts
    
    name = subs{subby} ;
    stimtrigcounts = [1,1,1] ;
    cd(['c:/shared/badger_eeg/',name,'/']) ; ls 
    setnames = dir('preproc*allstims*set') ; 
    eegsets = {setnames(1).name,setnames(2).name} ; 
    eeg = pop_loadset(eegsets{2}) ; 
    clear ersp ; 
    ep = pop_epoch(eeg,{'S  1','S  2','S  3'},[-2,7]) ; 
    for comp=1:size(ep.icaact,1) ; 
            [ersp(comp,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(comp,:,:)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                'plotitc','off','plotersp','off','baseline',0,'freqs',[1,128],'nfreqs',64,'winsize',64 );         
    end
   % comps = allcomps{subby} ; 
   figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(i,:,:)),[-4,4]) ; title(i) ; end


    cd(['C:\shared\badger_mri\',name,'\nii']) ; 
    bold = load_untouch_nii(fmrivols{2}) ; 
    boldim = bold.img ; 
    lowpower = log(eegfiltfft(eeg.icaact,eeg.srate,10,25).^2) ; 
    highpower = log(eegfiltfft(eeg.icaact,eeg.srate,35,70).^2) ; 
    lowpower = imresize(lowpower,[64,735]) ; 
    highpower = imresize(highpower,[64,735]) ; 
    wsize = floor(size(eeg.data,2)./size(boldim,4)) ; 
    clear stdpower 
    for c=1:64 ; icount = 1 ; 
        for i=1:wsize:size(eeg.data,2)-wsize
           stdpower(c,icount) = std(eeg.icaact(c,i:i+wsize)) ;
           icount = icount + 1 ; 
        end
    end
    
    clear convlow convhigh convstd
    for i=1:64 ;
        clow = conv(lowpower(i,:),hrf,'full') ; 
        convlow(i,:) = clow(1:size(lowpower,2)) ; 
        convlow(i,1:10) = squeeze(mean(convlow(i,:),2)) ; 
        chigh = conv(highpower(i,:),hrf,'full') ; 
        convhigh(i,:) = chigh(1:size(highpower,2)) ; 
        convhigh(i,1:10) = squeeze(mean(convhigh(i,:),2)) ; 
        cstd = conv(stdpower(i,:),hrf,'full') ; 
        convstd(i,:) = cstd(1:size(stdpower,2)) ; 
        convstd(i,1:10) = squeeze(mean(convstd(i,:),2)) ; 
    end
    
    
    
    f = load_untouch_nii('f_retino_allstims_01.nii') ; 
    for i=1:64 ; disp(i) ; 
        corrs(:,:,:,i) = voxcorr(boldim(:,:,:,1:700),squeeze(convlow(i,1:700))) ;       
    end
    zinds = 25:35; figure,
    for i=1:size(corrs,4) ; subplot(5,13,i),
        plotoverlayIntensity2D(squeeze(mean(f.img(zinds,:,1:end-2),1)),mat2gray(abs(squeeze(mean(corrs(zinds,:,1:end-2,i),1)))),squeeze(mean(corrs(zinds,:,1:end-2,i),1)),90) ;  title(i) ; 
    end

end
f = load_untouch_nii('f_retino_allstims_01.nii') ; 
mkdir corrs ; cd corrs ; 
for i=1:size(corrs,4);
    f.img = squeeze(corrs(:,:,:,i)) ; 
    save_untouch_nii(f,['corrs_',num2str(i),'.nii.gz']) ; 
end





