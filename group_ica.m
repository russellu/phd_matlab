clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','valerie'} ;
for subby=1:length(subs) ; 
    clear merged 
    cd(['c:/shared/badger_eeg/',subs{subby}]) ; %ls 
    prefix = 'setp_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ;gammas2=dir([prefix,'preproc*gamma*set']) ;
    allstims1=dir([prefix,'preproc*allstim*set']) ;allstims2=dir([prefix,'preproc*allstim*set']) ;
    movies=dir([prefix,'preproc*movie*set']) ; rests=dir([prefix,'preproc*rest*set']) ;
    %{
    if size(movies,1) ~= 0  
        setNames = {allstims1(1).name,allstims2(2).name,gammas1(1).name,gammas2(2).name,movies(1).name,rests(1).name} ; 
    else
        setNames = {allstims1(1).name,allstims2(2).name,gammas1(1).name,gammas2(2).name,rests(1).name} ; 
    end   
    %}
    setNames = {rests(1).name} ; 
    for sn=1:length(setNames)
       EEG = pop_loadset(setNames{sn}) ;  
       if sn==1 
           merged = EEG ;  
       else
           merged = pop_mergeset(merged,EEG) ;
       end
    end
    allmerge{subby} = merged ; 
end

for i=1:length(allmerge)
   if i==1 ; merged = allmerge{i} ; 
   else merged = pop_mergeset(merged,allmerge{i}) ; 
   end
end


mergefilt = merged ;
mergefilt.data = eegfiltfft(merged.data,merged.srate,1,128) ; 
mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,merged.srate,59.6,60.5) ; 

icamerge = pop_runica(mergefilt,'runica') ; 
pop_saveset(icamerge,'icamerge_rest') ; 
copy = icamerge ; 







