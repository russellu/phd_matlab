cd('C:\shared\badger\alex') ; 
EEG = pop_loadbv('.','alex_outside.vhdr') ; 
EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
EEG = pop_resample(EEG,256) ; 
filt = EEG ; 
filt.data = eegfiltfft(EEG.data,EEG.srate,1,128) ; 
filt = pop_runica(filt,'runica') ; 

filtepoch = pop_epoch(filt,{'S  1'},[-2,4]) ; 

ep = squeeze(filtepoch.icaact(26,:,3)) ; 
comp = squeeze(filt.icaact(26,:)) ; 
for i=1:length(comp)-length(ep)
    corrs(i) = corr2(ep,comp(i:i+length(ep)-1)) ; 
end


