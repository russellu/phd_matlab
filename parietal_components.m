clear all ; close all; 
cd E:\eegs\badger_alex

sets = dir('den*set'); 
for st=1:length(sets)
    eeg = pop_loadset('den_retino_rest.set'); 
    if st==1; merged= eeg; else merged = pop_mergeset(eeg,merged); end
end
filteeg = merged; filteeg.data = eegfiltfft(merged.data,merged.srate,1,128); 
[weights,sphere] = runica(filteeg.data,'maxsteps',128); 
winv = pinv(weights*sphere); 

for i=1:64; subplot(5,13,i) ; topoplot(winv(:,i),eeg.chanlocs); end