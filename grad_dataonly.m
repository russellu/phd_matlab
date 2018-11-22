clear all ; close all 
cd c:/shared/badger_eeg/alex ; 
mov = dir('*vhdr') ; 
%cd c:/users/acer/downloads/forRussell
%mov = dir('*edf') ;
for mm=1:length(mov)
EEG = pop_loadbv('.',mov(mm).name) ; 
%EEG = pop_readedf(mov(mm).name) ; 
%cd C:\Users\Acer\Downloads\ForRussell
%edf = dir('*edf') ; 
%EEG = pop_readedf(edf(4).name) ; 

EEGdata = EEG.data ; srate = 5000 ; 
% remove gradient from EEG dataset
chanN = 2 ;     
% step 1: get the indices of the gradient epochs 
hp = eegfiltfft(EEGdata(chanN,:),5000,20,5000) ; 
allhp = eegfiltfft(EEGdata([1:round(size(EEGdata,1)/2)],:),srate,30,round(srate/2)) ; 
stdhp = std(allhp,0,1) ; 
dhp = diff(hp,2,2) ; 
xc = xcorr(hp,srate*5,'coeff') ; 
wsize = median(diff(find(xc>.98))) ; % slice or TR window size (in samples) 
allsi = zeros(round(length(hp)./wsize),wsize) ; testepochs = zeros(size(allsi)) ; 
icount = 1 ; maxi = 100 ; 
for i=1:wsize:length(hp)-wsize*2 ;    
    allmaxes(icount) = maxi+i-1 ; 
    allsi(icount,:) = i+maxi:i+maxi+wsize-1 ; 
    testepochs(icount,:) = hp(allsi(icount,:)) ; 
    icount = icount + 1 ; 
end
zmaxes = zeros(1,size(EEG.data,2)) ; zmaxes(allmaxes) = 1 ; 


% step 2 find the indices where no gradient is present and remove those (ie, start and end effects):
sumsv = sum(diff(testepochs,1,2).^2,2) ; 
medsv = median(sumsv) ; stdsv = std(sumsv) ; 
bads = (sumsv < medsv - stdsv | sumsv > medsv + stdsv) ; 
allsi(bads,:) = [] ;  testepochs(bads,:) = [] ; teps = testepochs(:,1:10:end) ; 
corrs = corr(teps') ; bad2s = find(zscore(sum(corrs))<-6) ; 
allsi(bad2s,:) = [] ; testepochs(bad2s,:) = [] ; 
timep(1) = allsi(1,end) ; timep(2) = allsi(end,1) ; 

padsz = round(0.1*size(testepochs,1)) ; if padsz < 25 ; padsz = size(testepochs,1)-1 ; end % bound checking
if mod(padsz,2) == 1 ; padsz = padsz - 1 ; end % bound checking
if size(testepochs,1)-padsz <= 1 ; padsz = padsz - 2 ; end % bound checking

% step 3 subtract the average artifact from each channel separately: 
cleanchans = zeros(size(EEG.data)) ; 
for chanN=1:size(EEGdata,1) ; disp(['chan=',num2str(chanN)]) ; 
    chan = EEGdata(chanN,:) ; 
    fftchan = fft(chan) ; freqbinsz = srate/length(fftchan) ; 
    high_cutoff = 12 ; 
    lowp_ind = round(high_cutoff/freqbinsz) ; lpfft = fftchan ; lpfft(lowp_ind:end-lowp_ind) = 0 ; lowpower = real(ifft(lpfft)) ; 
    hpfft = fftchan ; hpfft([1:lowp_ind,end-lowp_ind:end]) = 0 ; highpower = real(ifft(hpfft)) ; 
    % for all gradient indices, get the artifact at that index:
    allhp = zeros(size(allsi)) ; allraw = zeros(size(allsi)) ; 
    for i=1:size(allsi,1)
        allraw(i,:) = highpower(allsi(i,:)) ; 
    end
    filt = fspecial('gaussian',[300,1],50) ; 
    padchan = zeros(size(allraw,1)+padsz*2,size(allraw,2)) ; padchan(padsz+1:end-padsz,:) = allraw ; padchan(1:padsz,:) = (flipud(allraw(1:padsz,:))) ; 
    padchan(end-padsz-1:end,:) = (flipud(allraw(end-padsz-1:end,:))) ; 
    filtchan = imfilter(padchan,filt,'replicate') ; filtchan = filtchan(padsz-1:end-padsz,:) ;
    subts = zeros(size(chan)) ; 
    for i=1:size(allsi,1) ; subts(allsi(i,:)) = allraw(i,:) - filtchan(i,:) ; end
    subts_final = subts + lowpower ; 
    cleanchans(chanN,:) = subts_final ; 
end 
clear reschans gradchans
% resample the final result
for i=1:size(cleanchans,1)
    reschans(i,:) = resample(cleanchans(i,:),250,srate) ; 
    gradchans(i,:) = resample(double(EEGdata(i,:)),250,srate) ; 
end
[s,f] = spectopo(reschans(:,15000:end-15000),0,250,'plot','off') ; 
[sg,fg] = spectopo(gradchans,0,250,'plot','off') ; 
figure,
subplot(2,2,1) ; imagesc(f,1:64,s) ; title('Removed (without gradient)') ; 
subplot(2,2,2) ; imagesc(f,1:64,sg) ; title('Original (with gradient') ;
subplot(2,2,3) ; imagesc(f,1:64,sg-s) ; title('Original minus Removed') ; xlabel('frequency(hz)') ; ylabel('channel #') ; 
subplot(2,2,4) ; plot(sg(5,:)) ; hold on ; plot(s(5,:),'r') ;title('Single channel, red=removed,blue=original') ; 
ylabel('power(db)') ; set(gca,'XTick',1:50:length(f),'XTickLabel',round(f(1:50:end))) ; xlabel('frequency(hz)') ; 
suptitle('Gradient denoising (All channels)') ; 
end

filtchans = eegfiltfft(reschans,250,1,250) ;
[weights,sphere] = runica(filtchans,'maxsteps',128) ;
acts = weights*sphere*filtchans ; 


