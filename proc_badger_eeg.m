cd c:/shared/badger/alex ; ls ; clear all ; close all ; 
vhdrs = dir('*retino*vhdr') ; 
for i=1:length(vhdrs)
    EEG = pop_loadbv('.',vhdrs(i).name) ; 
    EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
    EEG = denoise_grad2(EEG) ; 
    EEG = denoise_bcg(EEG) ; 
    if i==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end
    eegs{i} = EEG ; 
end

mergefilt = merged ; mergefilt.data = eegfiltfft(merged.data,merged.srate,1,128) ; 
mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,mergefilt.srate,59.5,60.5) ; 
ica = pop_runica(mergefilt,'runica') ; 

for i=1:length(eegs)
   eeg2 = eegs{i} ; 
   eeg2.icaact = icaact(eeg2.data,ica.icaweights*ica.icasphere,0) ; 
   eeg2.icachansind = ica.icachansind ; 
   eeg2.icasphere = ica.icasphere ; 
   eeg2.icasplinefile = ica.icasplinefile ; 
   eeg2.icaweights = ica.icaweights ; 
   eeg2.icawinv = ica.icawinv ;   
   eeg2 = pop_saveset(eeg2,[strrep(vhdrs(i).name,'.vhdr','')]) ; 
end



[s,f] = spectopo(eeg2.icaact,0,eeg2.srate,'plot','off') ; 
[s2,f2] = spectopo(eeg2.data,0,eeg2.srate,'plot','off') ; 
subplot(1,2,1) ; imagesc(f,1:64,s) ; title('component power') ; ylabel('component #') ; subplot(1,2,2) ; imagesc(f,1:64,s2) ; title('electrode power') ; ylabel('electrode #') ; xlabel('frequency(hz)') ; 











