clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
for sub=9%:length(subs)
cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls
mov = dir('bcgica*rest*set') ; 
EEG = pop_loadset(mov(1).name) ; 
EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
ica = pop_runica(EEG,'runica') ; 
end

freqs = 1:80 ;
allpows = zeros(size(ica.icaact,1),length(freqs),size(ica.icaact,2)) ; 
for f=1:length(freqs)
    allpows(:,f,:) = eegfiltfft(ica.icaact,250,freqs(f)-2,freqs(f)+2) ;   
end


evts = {EEG.urevent.type} ; 
lats = cell2mat({EEG.urevent.latency}) ; 
s1s = find(strcmp('S  1',evts)) ; 
r128s = find(strcmp('R128',evts)) ; 
s1lats = lats(s1s) ; r128lats = lats(r128s) ; 
ntrs = length(r128s)-1 ; 
cd(['c:/shared/newbadger_mri/',subs{sub}]) ; ls
bp = load_untouch_nii('bp_reg_topup_mc_retino_rest.nii.gz') ; bpimg = bp.img ; 

corr = load_untouch_nii('corrs.nii.gz') ; 
mask = medfilt3(corr.img > .12) ; inds = find(mask==1) ; [i1,i2,i3] = ind2sub(size(mask),inds) ; 
for i=1:length(inds) ; ts(i,:) = squeeze(bpimg(i1(i),i2(i),i3(i),:)) ; end ; mts = mean(ts,1) ; mts = mts(1:ntrs) ; 
%{
clear xcorrs 
for mc=1:64 ;
    mchans = squeeze(mean(allpows(mc,:,:).^2,1)) ; 
    mchans = eegfiltfft(mchans,250,0.02,2) ; 
    mchans = imfilter(mchans,fspecial('gaussian',[1,200],30)) ; 
    st = r128lats(1) ; en = r128lats(end) ; 
    newmd = mchans(:,st:en) ; resmd = imresize(newmd,[80,ntrs]) ; 
    for i=1:size(resmd,1)
       xcorrs(mc,i,:) = xcorr(mts(35:end-35),smooth(resmd(i,35:end-35)),20,'coeff') ; 
    end
end
for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(xcorrs(i,:,:)),[-1,1]) ; axis xy ; title(i) ; end

goods = [3,5,6,13,15,17,18,19,20,26,31] ; 
figure,imagesc(squeeze(mean(xcorrs(goods,:,:))),[-.4,.4]) ; axis xy ; 
%}


clear trchans
for mc=1:64 ;
    mchans = squeeze(mean(allpows(mc,:,:).^2,1)) ; 
    mchans = eegfiltfft(mchans,250,0.02,2) ; 
    mchans = imfilter(mchans,fspecial('gaussian',[1,200],30)) ; 
    st = r128lats(1) ; en = r128lats(end) ; 
    newmd = mchans(:,st:en) ; resmd = imresize(newmd,[80,ntrs]) ; 
    trchans(mc,:,:) = resmd ; 
end

restr = reshape(trchans,[64,numel(trchans(1,:,:))]) ; 
[U,S,V] = svd(restr') ; 
figure,
for i=1:50
res2 = reshape(U(i,:),[80,449]) ;
subplot(5,10,i) ; plot(log(mean((abs(res2'))))) ; 
end

clear xcorrs 
for mc=1:64 ;
    resmd = squeeze(resacts(mc,:,:)) ; 
    for i=1:size(resmd,1)
       xcorrs(mc,i,:) = xcorr(mts(35:end-35),smooth(resmd(i,35:end-35)),20,'coeff') ; 
    end
end
for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(xcorrs(i,:,:)),[-1,1]) ; axis xy ; title(i) ; end





