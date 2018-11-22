clear all ; close all ; 
cd c:/shared/resmerged/russell ; 
EEG = pop_loadset('merged.set') ; 
filt = eegfiltfft(EEG.data,EEG.srate,1,128) - eegfiltfft(EEG.data,EEG.srate,59,61) ; 
EEG.data = filt ; 
ica = pop_runica(EEG,'runica') ; 
acts = ica.icaweights*ica.icasphere*EEG.data ; 

subplot(2,2,1) ; plot(acts(2,50000:55000)) ; xlim([0,5000]) ; xlabel('sample (256Hz)')  ; ylabel('a.u.') ;
subplot(2,2,2) ; plot(acts(6,70000:71090)) ; xlim([0,1000]) ; ylim([-6,6]) ; xlabel('sample (256Hz)')  ; ylabel('a.u.') ;
subplot(1,2,1) ; topoplot(ica.icawinv(:,2),EEG.chanlocs) ; subplot(1,2,2) ; topoplot(ica.icawinv(:,6),EEG.chanlocs) ; 