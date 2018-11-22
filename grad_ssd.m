clear all ; close all ; 
%subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'} ; 
subs = {'MONG_01_RB','MONG_02_DP','MONG_03_CG','MONG_04_AB','MONG_05_SG','MONG_06_TS'} ; 
for subby=1%:length(subs)
cd(['c:/shared/mong_eeg/',subs{subby}]) ; ls 
name = dir(['*vhdr']) ; 
for nm=1%:length(name)
eeg = pop_loadbv('.',name(nm).name) ; 
eeg = pop_chanedit(eeg,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;

higheeg = eeg ; 
f = fft(higheeg.data,[],2) ;  
freqincr = 5000/(size(f,2)) ;
freqs = freqincr:freqincr:2500 ; 
freqaxis = zeros(1,size(f,2)) ; freqaxis(1:length(freqs)) = freqs ;
freqaxis(end-length(freqs)+1:end) = fliplr(freqs) ; 
raw_lowf = f(:,freqaxis<5) ; newf = zeros(size(f)) ; newf(:,freqaxis<5) = raw_lowf ; lowdat = real(ifft(newf,[],2)) ; 
f(:,freqaxis<5) = 0 ; higheeg.data = real(ifft(f,[],2)) ;
hfeeg = eegfiltfft(eeg.data,eeg.srate,40,2500) ; 

slice = 375 ; 
epochs = zeros(64,round(size(eeg.data,2)/slice),slice) ; rawepochs = zeros(64,round(size(eeg.data,2)/slice),slice) ; 
epochinds = zeros(round(size(eeg.data,2)/slice),slice) ; icount = 1 ; 
for i=1:slice:size(eeg.data,2)-slice
    rawepochs(:,icount,:) = higheeg.data(:,i:i+slice-1) ; 
    epochs(:,icount,:) = hfeeg(:,i:i+slice-1) ; 
    epochinds(icount,:) = i:i+slice-1 ; 
    icount = icount + 1 ; 
end

sumepochs = squeeze(mean(sum(epochs.^2,2),3)) ; 
%[~,csi] = sort(sumepochs,'descend') ; 
csi = [11,41,35,64,38,6,60,26,15] ; 
ssds = zeros(9,size(epochs,2),size(epochs,2)) ; 
for i=1:9 ; disp(i) ; 
    for j=1:size(epochs,2) ; disp(j) ; 
        ssds(i,j,:) =  sum(abs(repmat(squeeze(epochs(csi(i),j,1:5:end))',[size(epochs,2),1]) - squeeze(epochs(csi(i),:,1:5:end))),2) ;
    end
end
mssds = squeeze(mean(ssds(:,:,:),1)) ; clear ssds ; 
[sv,si] = sort(mssds,2,'ascend') ; 
% create the average epoch:
imagesc(si) ; 

newdat = zeros(size(eeg.data)) ; 
for i=1:size(mssds,1)
    newdat(:,epochinds(i,:)) = squeeze(rawepochs(:,i,:)) - squeeze(mean(rawepochs(:,si(i,2:40),:),2)) ; 
end

postepochs = zeros(64,round(size(eeg.data,2)/slice),slice) ; icount = 1 ; 
for i=1:slice:size(eeg.data,2)-slice
    postepochs(:,icount,:) = newdat(:,i:i+slice-1) ; 
    icount = icount + 1 ; 
end

fulldat = newdat + lowdat ; 

end
end