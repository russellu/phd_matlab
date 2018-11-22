cd('c:/shared/sim_1/first_simultaneous_recordings') ; ls 
clear all ; close all ; 
for eset=1:5
EEG = pop_loadbv('.',['EEG-fMRI_Russell_visual_',num2str(eset),'.vhdr']) ;
%EEG = pop_loadbv('.',['EEG-fMRI_Russell_',num2str(eset),'.vhdr']) ; 
EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
mrindices = find(strcmp({EEG.urevent.type},'R128')) ; 
lats = cell2mat({EEG.urevent.latency}) ; gradlats = lats(mrindices) ; 
%topoplot(1:64,EEGdenoise.chanlocs,'electrodes','numbers')
%c1  = [46,52,60] ; c2 = [11,47,53] ; 
EEG2 = EEG ; 
for chan=1:64 ;
diffgradlats = diff(gradlats) ;
gradepochs = zeros(length(gradlats),diffgradlats(2)) ; 
for i=1:length(gradlats)
        if gradlats(i)+diffgradlats(1) < size(EEG.data,2) ; 
        gradepochs(i,:) = EEG.data(chan,gradlats(i):gradlats(i)+diffgradlats(1)-1) ; 
        end
    
end
mepochs = squeeze(mean(gradepochs)) ; 
for i=1:length(gradlats)-1
     EEG2.data(chan,gradlats(i):gradlats(i)+diffgradlats(2)-1) =  EEG.data(chan,gradlats(i):gradlats(i)+diffgradlats(2)-1)-mepochs ; 
end
end
res = pop_resample(EEG2,500) ; 


%[t,f] = spectopo(res.data(46,20000:180000),0,res.srate,'plot','off') ; 
EEG3 = pop_select(EEG2,'point',[gradlats(2),gradlats(length(gradlats)-1)]) ; 
res = pop_resample(EEG3,500) ; 

%EEG3 = pop_resample(EEG3,500) ; 
if eset==1
    merged = res ; 
else merged = pop_mergeset(res,merged) ; 
end
end

%merged2 = merged ; merged2.data = eegfiltfft(merged2.data,merged2.srate,1,120) ; 
%merged2 = pop_runica(merged2,'runica') ; 
%eeg = pop_loadset('merged_russell_visual.set') ; 
%both  = pop_mergeset(eeg,merged) ; 
%pop_saveset(both,'merged_russell_dave') ; 



merged = pop_saveset(merged,'merged_russ_1_5.set') ; 
merged2 = pop_loadset('merged_dave_1_5.set') ; 
both = pop_mergeset(merged,merged2) ; 
both = pop_saveset(both,'bothmerged') ; 



