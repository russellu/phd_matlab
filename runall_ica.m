subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'} ; 
for s=2:length(subs) ; 
cd c:/shared/badger_eeg2 ; ls ; cd(subs{s}) ; ls 
allsets = dir('bcgica*set') ;
for set=1:length(allsets) 
    EEG = pop_loadset(allsets(set).name) ; 
    if set==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG) ; end
end
ica = pop_runica(merged,'runica') ; 
ica = pop_chanedit(ica,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
allica{1} = ica.icaweights ; allica{2} = ica.icasphere ; 
save('allica','allica') ; 
end




subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
for s=1:length(subs) ; 
cd c:/shared/badger_eeg2 ; ls ; cd(subs{s}) ; ls 
allsets = dir('bcgica*set') ;
    EEG = pop_loadset(allsets(1).name) ; 
EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;

allica = load('allica') ; allica = allica.allica ; 
weights = allica{1} ; sphere = allica{2} ; 
winv = pinv(weights*sphere) ; 
figure ; for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),EEG.chanlocs) ; title(i); end ; suptitle(subs{s}) 
end

comps = {[5,12,14,16,21,19],
    [6,14,16,17,25],
    [6,10,111,13,16,28],
    [6,9,23],
    [12,15,18],
    [8,15,20],
    [7,13,21,31],
    [6,15,21,22],
    [8,15,18]} ;



