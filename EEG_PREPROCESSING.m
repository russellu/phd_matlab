%% PP3

nb_participants = 1;

raw_eeg_data = dir('*.vhdr');

for i = 1:nb_participants

subject_code = raw_eeg_data(i).name
    
EEG = pop_loadbv('R:\\Neuro\\Rech_PMBernier\\EtudeEEG-IRM\\Felix\\Laphroaig\\EEG_DATA', subject_code, [], [1:64]);

EEG=pop_chanedit(EEG, 'lookup','C:\\Program Files\\MATLAB\\R2014a\\eeglab12_0_2_2b\\plugins\\dipfit2.2\\standard_BESA\\standard-10-5-cap385.elp');

EEG = pop_resample(EEG, 256);

%filter... Ask Russ and KW...

%Bandpass filter (usage) -> [smoothdata] = eegfiltfft(data,srate,locutoff,hicutoff);
EEG.data = eegfiltfft(EEG.data, EEG.srate, 0.1, 100);

EEG = pop_reref( EEG, []);

%EEG = pop_saveset( EEG, 'filename',[int2str(i) '.set'],'filepath','R:\Neuro\Rech_PMBernier\EtudeEEG-IRM\Felix\Laphroaig\EEG_DATA\EEG_DATA_PP3');

end
