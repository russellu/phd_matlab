clear all ; close all ; 
cd C:\shared\muthu_rest

EEG = pop_loadbv('.','EFKAM008_PLA.vhdr') ; 
evts = {EEG.urevent.type} ; 
lats = cell2mat({EEG.urevent.latency}) ; 
trigs = {'V  1'} ; 
gradinds = find(strcmp('V  1',evts)) ; 

gradlats = lats(gradinds) ; 
difflats = diff(gradlats) ; 
modedifflats = mode(difflats) ; 
transitions = find(difflats ~= modedifflats) + 1 ; 
newtrans = zeros(1,length(transitions)+1) ; newtrans(2:end) = transitions ; newtrans(1) = 1 ; newtrans(end+1) = (length(gradlats)) ; 

for i=1:length(newtrans)-1
    onsetinds{i} = newtrans(i):newtrans(i+1)-1 ; 
end


indsi = onsetinds{5} ; 
latsi = gradlats(indsi) ; 
st = latsi(1)-40000 ; en = latsi(end)+40000 ; 

eeg2 = pop_select(EEG,'point',[st,en]) ; 


EEG = eeg2 ; 

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
    padchan = zeros(size(allraw,1)+100,size(allraw,2)) ; padchan(51:end-50,:) = allraw ; padchan(1:50,:) = (flipud(allraw(1:50,:))) ; 
    padchan(end-49:end,:) = (flipud(allraw(end-49:end,:))) ; 
    filtchan = imfilter(padchan,filt,'replicate') ; filtchan = filtchan(51:end-50,:) ; 
    subts = zeros(size(chan)) ; 
    for i=1:size(allsi,1) ; subts(allsi(i,:)) = allraw(i,:) - filtchan(i,:) ; end
    subts_final = subts + lowpower ; 
    cleanchans(chanN,:) = subts_final ; 
end 
EEG2 = EEG ; 
EEG2.data = double(cleanchans) ; 
EEGres = pop_resample(EEG2,250) ;   timep = timep./20 ; 
clear reschans 
for i=1:size(cleanchans,1)
   reschans(i,:) = resample(cleanchans(i,:),250,5000) ;  
end
reschans(32,:) = rand(1,length(reschans))*0.0005 ; 
[s,f] = spectopo(reschans,0,250,'plot','off') ; 
filtreschans = eegfiltfft(reschans,250,1,250) ; 
[weights,sphere] = runica(filtreschans,'maxsteps',100) ; 
acts = weights*sphere*filtreschans ; 
[sica,fica] = spectopo(acts,0,250,'plot','off') ; 
EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;

[newz,eegbcg] = fbcg(EEGres) ; 
[weights,sphere] = runica(newz,'maxsteps',100) ; 
winv = pinv(weights*sphere) ; 
[ss,ff] = spectopo(weights*sphere*newz,0,250,'plot','off') ; 



