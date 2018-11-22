clear all ; close all ; 
cd C:\shared\discret
discr = dir('grating*vhdr') ; 
for i=1:length(discr)
   EEG = pop_loadbv('.',discr(i).name) ; 
   EEG = pop_resample(EEG,256) ; 
   if i==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG) ; end
end
merged = pop_chanedit(merged,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
temp = merged.data(1:32,:) ; merged.data(1:32,:) = merged.data(33:64,:) ; merged.data(33:64,:) = temp ; 
mergefilt = merged ; mergefilt.data = eegfiltfft(merged.data,merged.srate,40,80) - eegfiltfft(merged.data,merged.srate,59,61) - eegfiltfft(merged.data,merged.srate,119,121) ; 
ica = pop_runica(mergefilt,'runica') ; ica.icaact = ica.icaweights*ica.icasphere*merged.data ; 
[s,f] = spectopo(ica.icaact,0,ica.srate,'plot','off') ; 


