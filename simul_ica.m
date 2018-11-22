cd c:/shared/simul_EEG ; ls 
EEG = pop_loadbv('.','stim_01.vhdr')  ;
EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
EEG.data = eegfiltfft(EEG.data,EEG.srate,100,2500) ; 
EEG = pop_runica(EEG,'runica') ; 


