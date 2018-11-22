clear all ; close all ; 
%subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'} ; 
subs = {'MONG_01_RB','MONG_02_DP','MONG_03_CG','MONG_04_AB','MONG_05_SG','MONG_06_TS'} ; 
for subby=3%:length(subs)
cd(['c:/shared/mong_eeg/',subs{subby}]) ; ls 
name = dir(['*vhdr']) ; 
for nm=3%:length(name)
eeg = pop_loadbv('.',name(nm).name) ; 
eeg = pop_chanedit(eeg,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;

higheeg = eeg ; 
f = fft(higheeg.data,[],2) ;  
freqincr = 5000/(size(f,2));
freqs = freqincr:freqincr:2500 ; 
freqaxis = zeros(1,size(f,2)) ; freqaxis(1:length(freqs)) = freqs ;
freqaxis(end-length(freqs)+1:end) = fliplr(freqs) ; 
raw_lowf = f(:,freqaxis<5) ; newf = zeros(size(f)) ; newf(:,freqaxis<5) = raw_lowf ; lowdat = real(ifft(newf,[],2)) ; 
f(:,freqaxis<5) = 0 ; higheeg.data = real(ifft(f,[],2)) ;

slice = 375 ; epochs = zeros(64,round(size(eeg.data,2)/slice),slice) ; icount = 1 ; 
epochinds = zeros(round(size(eeg.data,2)/slice),slice) ; 
for i=1:slice:size(eeg.data,2)-slice
    epochs(:,icount,:) = higheeg.data(:,i:i+slice-1) ; 
    epochinds(icount,:) = i:i+slice-1 ; 
    icount = icount + 1 ; 
end

sumepochs = zscore(squeeze(mean(sum(abs(epochs),3),1))) ;
bads = find(sumepochs<-1) ; 
goods = find(sumepochs>-1) ; 
rndinds = randperm(length(goods)) ; rndinds = goods(rndinds(1:length(bads))) ; 
epochs(:,bads,:) = epochs(:,rndinds,:) ; 
alldiffs = zeros(size(epochs)) ; 
for i=1:size(epochs,1) ; disp(i) ;
    epochsi = squeeze(epochs(i,:,:)) ; 
    pad = zeros(size(epochsi,1)+200,size(epochsi,2)) ; 
    pad(101:end-100,:) = epochsi ; midind = round(size(pad,1)/2) ; 
    pad(1:100,:) = pad(midind:midind+99,:) ; pad(end-100:end,:) = pad(midind:midind+100,:) ;
    smthepochs = imfilter(pad,fspecial('gaussian',[160,1],60)) ; 
    diffs = pad-smthepochs ; 
    diffs = diffs(101:end-100,:) ; 
    alldiffs(i,:,:) = diffs ; 
end
newts = higheeg.data ; 
epochinds(epochinds==0) = 1 ; 
for i=1:size(alldiffs,2)
    newts(:,epochinds(i,:)) = squeeze(alldiffs(:,i,:)) ; 
end
newts = newts + lowdat  ;
fnew = abs(fft(newts(25,200000:end-200000))) ; 
plot(fnew(3000:end/8)) ; hline(1.22*10^5) ; 

epochs = zeros(64,round(size(eeg.data,2)/slice),slice) ; icount = 1 ; 
for i=1:slice:size(eeg.data,2)-slice
    epochs(:,icount,:) = newts(:,i:i+slice-1) ; 
    icount = icount + 1 ; 
end

end
end