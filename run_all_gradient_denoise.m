clear all ; close all ; 
% problem subjects: sukhman, tegan, jeremie, genevieve. 
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'} ;
subcount = 1 ; 
for s=1:length(subs) ; 
    cd(['c:/shared/badger_eeg/',subs{s}]) ; ls *vhdr ;
    vhdrs = dir('*vhdr') ; 
    for vhdr=1:length(vhdrs)
        EEG = pop_loadbv('.',vhdrs(vhdr).name) ; % load the data 
        [EEGnograd,timep] = remove_gradient2(EEG) ; timep = round(timep) ;
        [spec(subcount,vhdr,:,:),f] = spectopo(EEGnograd.data(:,5000:end-5000),0,250,'plot','off') ; 
        %EEGnograd = pop_chanedit(EEGnograd,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
        pop_saveset(EEGnograd,['subgrad_',strrep(vhdrs(vhdr).name,'.vhdr','')]) ; 
    end
    subcount = subcount + 1 ; 
end


