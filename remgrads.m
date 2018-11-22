clear all ; close all; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'}; 
for sb=1:length(subs)
cd(['E:\badger_eeg\',subs{sb}]) ; 
alleegs = dir('*vhdr');
for vhdr=1:length(alleegs)

eeg = pop_loadbv('.',alleegs(vhdr).name); 
eeg = pop_chanedit(eeg,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
newsrate = 5000 ; 
higheeg = eeg ; 
f = fft(higheeg.data,[],2) ;  
freqincr = 5000/(size(f,2));
freqs = freqincr:freqincr:2500 ; 
freqaxis = zeros(1,size(f,2)) ; freqaxis(1:length(freqs)) = freqs ;
freqaxis(end-length(freqs)+1:end) = fliplr(freqs) ; 
raw_lowf = f(:,freqaxis<5) ; 
lowfreqs = f; lowfreqs(:,freqaxis>=5) = 0 ; lowfreqs = real(ifft(lowfreqs,[],2)); 
f(:,freqaxis<5) = 0 ; higheeg.data = real(ifft(f,[],2)) ;


lats = cell2mat({eeg.urevent.latency}); 
types = {eeg.urevent.type}; 
r128s = find(strcmpi('R128',types)); 
startind = lats(r128s(1))-315*11;

timeinds = startind:lats(end); 
epochs = zeros(64,length(r128s)*11-1,315); 
epochinds = zeros(length(r128s)*11-1,315); 
icount = 1 ;
for i=startind:315:startind+315*11*length(r128s)-1
   epochs(:,icount,:) = higheeg.data(:,i:i+314); 
   epochinds(icount,:) = i:i+314; 
   icount = icount + 1; 
end

padepochs = zeros(size(epochs,1),size(epochs,2)+100,size(epochs,3)); 
padepochs(:,51:end-50,:) = epochs; 
padepochs(:,1:50,:) = epochs(:,50:-1:1,:); 
padepochs(:,end-50:end,:) = epochs(:,end:-1:end-50,:); 
subepochs = zeros(size(padepochs)); 
for i=1:64
    imgi = squeeze(padepochs(i,:,:));
    smoothi = imfilter(imgi,fspecial('gaussian',[90,1],150)); 
    subepochs(i,:,:) = squeeze(padepochs(i,:,:)) - smoothi;  
end
subepochs = subepochs(:,51:end-50,:); 

for i=1:size(epochinds,1)
    higheeg.data(:,epochinds(i,:)) = subepochs(:,i,:); 
end

higheeg.data = higheeg.data + lowfreqs; 
higheeg.data(:,[1:startind,lats(r128s(end)):end]) = 0; 
higheeg = pop_resample(higheeg,250);
pop_saveset(higheeg,strrep(alleegs(vhdr).name,'.vhdr','.set')); 
plot(higheeg.data(12,:)); 
end
end






