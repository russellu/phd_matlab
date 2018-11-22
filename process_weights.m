clear all ; close all ; 
icamerge = pop_loadset('c:/shared/badger_eeg/valerie/icamerge_rest.set')  ;
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','valerie'} ;
for subby=1:length(subs) ; 
    clear merged 
    cd(['c:/shared/badger_eeg/',subs{subby}]) ; %ls 
    preprocs = dir('preproc*set') ; 
    for p=1:length(preprocs)
        setp = pop_loadset(preprocs(p).name) ;  
        setp.icact = icaact(setp.data,icamerge.icaweights*icamerge.icasphere,0) ; 
        setp.icachansind = icamerge.icachansind ; 
        setp.icasphere = icamerge.icasphere ; 
        setp.icasplinefile = icamerge.icasplinefile ; 
        setp.icaweights = icamerge.icaweights ; 
        setp.icawinv = icamerge.icawinv ; 
        setp.data = setp.data - eegfiltfft(setp.data,icamerge.srate,59.5,60.5) ; 
        setp.icaact = setp.icaact - eegfiltfft(setp.icaact,icamerge.srate,59.5,60.5) ; 
        pop_saveset(setp,['setrest_',preprocs(p).name]) ; 
    end   
end