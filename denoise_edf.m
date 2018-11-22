cd C:\shared\muthu_rest\ForRussell
EEG =  pop_readedf('s1_send_run01-edf.edf') ; 

chanN = 22 ;     
% step 1: get the indices of the gradient epochs 
hp = eegfiltfft(EEG.data(chanN,:),5000,20,5000) ; 
allhp = eegfiltfft(EEG.data([1:16],:),5000,30,5000) ; 
stdhp = std(allhp,0,1) ; 
dhp = diff(hp,2,2) ; 
xc = xcorr(hp,20000,'coeff') ; 
wsize = median(diff(find(xc>.98))) ; 
allsi = zeros(round(length(hp)./wsize),wsize) ; testepochs = zeros(size(allsi)) ; 

icount = 1 ; 
for i=1:wsize:length(hp)-wsize*2 ; 
    epochi = stdhp(i:i+wsize-1) ; % maxi = find(epochi==max(epochi),1) ; 
    maxi = 100 ; 
    allmaxes(icount) = maxi+i-1 ; 
    allsi(icount,:) = i+maxi:i+maxi+wsize-1 ; 
    testepochs(icount,:) = hp(allsi(icount,:)) ; 
    icount = icount + 1 ; 
end
zmaxes = zeros(1,size(EEG.data,2)) ; zmaxes(allmaxes) = 1 ; 


% step 2 find the indices where no gradient is present and remove those:
sumsv = sum(diff(testepochs,1,2).^2,2) ; 
medsv = median(sumsv) ; stdsv = std(sumsv) ; 
%plot(sumsv) ; hline(medsv,'k') ; hline(medsv+stdsv,'r') ; hline(medsv-stdsv,'r') ; 
bads = (sumsv < medsv - stdsv | sumsv > medsv + stdsv) ; 
allsi(bads,:) = [] ;  testepochs(bads,:) = [] ; teps = testepochs(:,1:10:end) ; 
corrs = corr(teps') ; bad2s = find(zscore(sum(corrs))<-6) ; 
allsi(bad2s,:) = [] ; testepochs(bad2s,:) = [] ; 
timep(1) = allsi(1,end) ; timep(2) = allsi(end,1) ; 

% step 3 subtract the average artifact from each channel separately: 
cleanchans = zeros(size(EEG.data)) ; 
for chanN=1:size(EEG.data,1) ; disp(['chan=',num2str(chanN)]) ; 
    chan = EEG.data(chanN,:) ; 
    fftchan = fft(chan) ; freqbinsz = 5000/length(fftchan) ; 
    high_cutoff = 5 ; 
    lowp_ind = round(high_cutoff/freqbinsz) ; lpfft = fftchan ; lpfft(lowp_ind:end-lowp_ind) = 0 ; lowpower = real(ifft(lpfft)) ; 
    hpfft = fftchan ; hpfft([1:lowp_ind,end-lowp_ind:end]) = 0 ; highpower = real(ifft(hpfft)) ; 
    % for all gradient indices, get the artifact at that index:
    allhp = zeros(size(allsi)) ; allraw = zeros(size(allsi)) ; 
    for i=1:size(allsi,1)
        allraw(i,:) = highpower(allsi(i,:)) ; 
    end
    filt = fspecial('gaussian',[80,1],80) ; 
    padchan = zeros(size(allraw,1)+200,size(allraw,2)) ; padchan(101:end-100,:) = allraw ; padchan(1:100,:) = (flipud(allraw(1:100,:))) ; 
    padchan(end-99:end,:) = (flipud(allraw(end-99:end,:))) ; 
    filtchan = imfilter(padchan,filt,'replicate') ; filtchan = filtchan(101:end-100,:) ; 
    subts = zeros(size(chan)) ; 
    for i=1:size(allsi,1) ; subts(allsi(i,:)) = allraw(i,:) - filtchan(i,:) ; end
    subts_final = subts + lowpower ; 
    cleanchans(chanN,:) = subts_final ; 
end 
EEG2 = EEG ; 
EEG2.data = double(cleanchans) ; 
EEGres = pop_resample(EEG2,250) ;   timep = timep./20 ; 

res = cleanchans(:,1:20:end) ; 
[s,f] = spectopo(res,0,250,'plot','off') ; 
for i=1:size(cleanchans,1)
   reschans(i,:) = resample(cleanchans(i,:),250,5000) ;  
end
[s,f] = spectopo(reschans,0,250,'plot','off') ; 
filtreschans = eegfiltfft(reschans,250,1,250) ; 
[weights,sphere] = runica(filtreschans) ; 
acts = weights*sphere*filtreschans ; 
[sica,fica] = spectopo(acts,0,250,'plot','off') ; 


