clear all ; close all ;
cd E:\gamma_deprivation

eeg = pop_loadbv('.','lyes_md_01.vhdr'); 
eeg = pop_resample(eeg,256) ;
filteeg = eegfiltfft(eeg.data,eeg.srate,1,100); 
[weights,sphere] = pop_runica(filteeg,'maxsteps',128); 