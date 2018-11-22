clear all ; close all; 
cd E:\badger_eeg\dina
eeg2 = pop_loadbv('.','retino_gamma_02.vhdr'); 
reseeg2 = pop_resample(eeg2,250); 

gradeeg = remove_gradient2(eeg2); 
gradeeg = pop_chanedit(gradeeg,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;

filt = eegfiltfft(gradeeg.data,gradeeg.srate,1,128); % high pass
bades = [61,29,62,30,55,11,53,1,63,2,54,12,56,14,58,16,60,10,20,9,59,15,57,13]; 

[weights,sphere] = runica(filt(bades,:),'maxsteps',128); 
acts = weights*sphere*filt(bades,:); 



