clear all ; close all ; 
cd C:\shared\simdenoise\sukhman
vhdrs = dir('*vhdr') ; 
for v=1:length(vhdrs) ; 
eeg = pop_loadbv('.',vhdrs(v).name) ;   
eeg = pop_chanedit(eeg,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
fdat = fft(eeg.data,[],2) ; 
cutoffhz = 50 ; srate = 5000 ; 
freqstep = srate/size(fdat,2) ; 
cutoffind = round(cutoffhz/freqstep) ; 
highfdat = fdat ; highfdat(:,[1:cutoffind,end-cutoffind:end]) = 0 ; 
invfdat = real(ifft(highfdat,[],2)) ; 

epl = 315 ; 
epochs = zeros(size(invfdat,1),ceil(size(invfdat,2)/epl),epl) ; 
icount = 1 ; 
for i=1:epl:size(invfdat,2)-epl
    epochs(:,icount,:) = invfdat(:,i:i+epl-1) ; 
    icount = icount + 1 ; 
end

zdiffs = (zscore(squeeze(mean(sum(abs(diff(epochs,1,3)),3),1)))) ;
%figure,plot(zdiffs) ; vline(find(zdiffs<-2)) ; 
%zdiffs(zdiffs<-1) = [] ; 
figure,plot(zdiffs) ; ylim([-2,2]) ;


end
% using changes in the gradient wave form to detect motion?









