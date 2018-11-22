function EEG = ica_applyweights(raw_eeg,filtered_eeg) 
% function EEG = ica_applyweights(raw_eeg,filtered_eeg) 
% apply the ICA weights from filtered_eeg to the raw_eeg dataset
EEG = raw_eeg ; 
EEG.icaact = filtered_eeg.icaweights*filtered_eeg.icasphere*raw_eeg.data ;  
EEG.icawinv = filtered_eeg.icawinv ; EEG.icaweights = filtered_eeg.icaweights ; EEG.icasphere = filtered_eeg.icasphere ;
EEG.icachansind = filtered_eeg.icachansind ; EEG.icasplinefile = filtered_eeg.icasplinefile ; 

end