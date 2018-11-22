clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
for sub=2:length(subs)
cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls
movs = dir('bcgica*set') ; 
for m=1:length(movs)
    EEG = pop_loadset(movs(m).name) ; 
    if m==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG) ; end
end

ica = pop_runica(merged,'runica') ; 
ica = pop_chanedit(ica,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp')  ; 

ica.data = [] ; ica.icaact = [] ; 
icaw{1} = ica.icaweights ; icaw{2} = ica.icasphere ; 
save('icaw','icaw') ; 

end