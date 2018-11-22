clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
for sub=6:9
    cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls
    rest = dir('bcgica*rest*set') ; 
    movie = dir('bcgica*movie*set') ; 
    rest = pop_loadset(rest(1).name) ; rest = pop_chanedit(rest,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
    movie = pop_loadset(movie(1).name) ; movie = pop_chanedit(movie,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
    merged = pop_mergeset(rest,movie) ; 
   
    
    ica = pop_runica(merged,'runica') ; 
    
    movie = ica_applyweights(movie,ica) ; 
    rest = ica_applyweights(rest,ica) ; 
    
    pop_saveset(movie,'weights_movie') ; 
    pop_saveset(rest,'weights_rest') ; 
end