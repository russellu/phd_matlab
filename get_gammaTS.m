clear all ; close all

subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie'} ; 
prefix='1hz_' ; 

for sb=1:length(subs) ; 
    cd(['c:/shared/badger_eeg/',subs{sb}]) ; ls 
    gammas=dir([prefix,'preproc*gamma*set']) ; 
    figure,
    for g=1:2 
        clear stims
        cd(['c:/shared/badger_eeg/',subs{sb}]) ; ls 
        EEG = pop_loadset(gammas(g).name) ; 
        events = {EEG.urevent.type} ; 
        lats = {EEG.urevent.latency} ;
        stims{1} = find(strcmp('S  1',events)) ; stims{2} = find(strcmp('S  2',events)) ; stims{3} = find(strcmp('S  3',events)) ; 
        allinds = [cell2mat(stims(1)),cell2mat(stims(2)),cell2mat(stims(3))] ; 
        [~,si] = sort(allinds) ; 
        stimlats = cell2mat(lats(allinds(si))) ; 

        fmris = find(strcmp('R128',events)) ; 
        fmrilats = cell2mat(lats(fmris)) ; 

        rawstimlats = stimlats - fmrilats(1) ; 

        % find the fmri TR indices to which each stimulus is closest:
        design = zeros(1,735) ; 
        longdesign = zeros(1,size(EEG.data,2)+round(EEG.srate/.693)) ; 
        for i=1:length(stimlats)

           longdesign(rawstimlats(i):rawstimlats(i)+5*EEG.srate) = 1 ; 
           lati = find(abs(fmrilats-stimlats(i))==min(abs(fmrilats-stimlats(i)))) ;  
           design(lati:lati+round(5/0.693)) = 1 ; 
        end
        shortdesign = imresize(longdesign,[1,735]) ; 
        cd(['c:/shared/badger_mri/',subs{sb},'/nii',]) ; ls 
        dlmwrite(['gammadesign',num2str(g),'.txt'],round(shortdesign(40:695))') ; 
        plot(shortdesign) ; hold on ; 
    end
    
    % write the allstims
    
    trigs = {'S 11','S 12','S 13','S 14','S 21','S 22','S 23','S 24','S 31','S 32','S 33','S 34','S 41','S 42','S 43','S 44','S 71','S 72','S 73','S 74','S 81','S 82','S 83','S 84'} ;

    cd(['c:/shared/badger_eeg/',subs{sb}]) ; ls 
    allstims=dir([prefix,'preproc*allstim*set']) ; 
    for g=1:2 
        clear stims
        cd(['c:/shared/badger_eeg/',subs{sb}]) ; ls 
        EEG = pop_loadset(allstims(g).name) ; 
        events = {EEG.urevent.type} ; 
        lats = {EEG.urevent.latency} ;
       
        allinds = find(ismember(events,trigs)) ; 
        
        [~,si] = sort(allinds) ; 
        stimlats = cell2mat(lats(allinds(si))) ; 

        fmris = find(strcmp('R128',events)) ; 
        fmrilats = cell2mat(lats(fmris)) ; 

        rawstimlats = stimlats - fmrilats(1) ; 

        inds = find(ismember(events,trigs{1})) ; 
        
        % find the fmri TR indices to which each stimulus is closest:
        design = zeros(1,735) ; 
        longdesign = zeros(1,size(EEG.data,2)+round(EEG.srate/.693)) ; 
        for i=1:length(stimlats)

           longdesign(rawstimlats(i):rawstimlats(i)+5*EEG.srate) = 1 ; 
           lati = find(abs(fmrilats-stimlats(i))==min(abs(fmrilats-stimlats(i)))) ;  
           design(lati:lati+round(5/0.693)) = 1 ; 
        end
        shortdesign = imresize(longdesign,[40,695]) ; 
        cd(['c:/shared/badger_mri/',subs{sb},'/nii',]) ; ls 
        dlmwrite(['allstimdesign_',num2str(g),'.txt'],shortdesign') ; 
        plot(shortdesign) ; hold on ; 
    end
  
end







