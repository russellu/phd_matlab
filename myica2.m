clear all ; close all ; 
cd c:/shared/kraken_Pilot/ ; ls 
EEG = pop_loadbv('.','M1.vhdr') ; 
EEG = pop_resample(EEG,250) ;
EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
filtEEG = EEG ; % copy the dataset, to preserve the unfiltered data for later use
filtEEG.data = filtEEG.data - eegfiltfft(filtEEG.data,filtEEG.srate,59.5,60.5) ; % filter out line noise
filtEEG.data = eegfiltfft(filtEEG.data,filtEEG.srate,1,125) ; % high pass filter above 1Hz
trigs{1} = {'S  2','R  2','S  8'} ;
trigs{2} = {'S 32','R  8','S128'} ; 
trigs{3} = {'S 10','R 32','S 34'} ; 
trigs{4} = {'S130','R128','S 40'} ; 
ep = pop_epoch(filtEEG,[trigs{1},trigs{2},trigs{3},trigs{4}],[-2,2]) ; % all the epochs at once, for the decomposition
ica = pop_runica(ep,'runica','extended',1) ;