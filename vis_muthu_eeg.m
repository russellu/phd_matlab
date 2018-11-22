clear all ; close all; 
cd E:\muthu_eeg
eeg = pop_loadbv('.','EFKAM008_PLA.vhdr'); 

subplot(2,1,1); 
plot(eeg.data(48,6340000:6340000+5000)); ylim([-14000,11000]); ylabel('\muV') ;xlabel('sample (HZ=5000)'); title('tr=2200ms, slice gap'); vline(1800,'r'); 
types = {eeg.urevent.type}; r128s = find(strcmpi('V  1',types)); lats = cell2mat({eeg.urevent.latency}); tr = median(diff(lats(r128s)))/5000; 


cd E:\badger_eeg\alex
eeg2 = pop_loadbv('.','retino_allstims_02.vhdr'); 
subplot(2,1,2); 
plot(eeg2.data(48,1000000:1000000+5000)); ylim([-5000,5000]); ylabel('\muV') ;xlabel('sample (HZ=5000)'); title('TR=693ms, no slice gap'); 
set(gcf,'Position',[100 100 500 400])
