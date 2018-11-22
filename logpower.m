clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 

comps = {[30,26],[32,23],[21,28],[47,46],[19,37],[47,38],[44,46]} ; % all subjects, right and left
fcomps = {[64,80,12],[41,31,40],[61,48,15],[6,84,79],[79,80,50],[75,77,41],[48,50]} ; 

for sub=1%:length(subs)

    %%%% EEG processing:
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    prefix = 'highfreq_' ;
    gammas1=dir([prefix,'preproc*gamma*set']) ;gammas2=dir([prefix,'preproc*gamma*set']) ;
    allstims1=dir([prefix,'preproc*allstim*set']) ;allstims2=dir([prefix,'preproc*allstim*set']) ;
    movies=dir([prefix,'preproc*movie*set']) ; rests=dir([prefix,'preproc*rest*set']) ;
    setNames = {allstims1(1).name,allstims1(2).name,gammas1(1).name,gammas1(2).name,movies(1).name,rests(1).name} ;   
    
    rest = pop_loadset(rests(1).name) ; 
    restsz = size(rest.icaact,2) ; 
    restacts = rest.icaact(:,restsz/8:restsz-restsz/8) ; 
    [s,f] = spectopo(restacts,0,rest.srate,'plot','off') ; 
    
end