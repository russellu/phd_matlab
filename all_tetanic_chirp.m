clear all ; close all ; 
stims = {'S  1','S  2','S 11','S 12'} ;
subs = {'AD','AM','CV','DT','FMP','KC','LB'};
bades = {[],[],[21,61],[9,42],[19,57,62],[1,32,33,37],[10]};

cd(['E:\jf_data\Tetanic Visual\RB']) ; rmerged = pop_loadset('merged.set'); 


for sb=3%:length(subs)
    
    cd(['E:\jf_data\ERP_Chirp\',subs{sb}]) ;
    
    %{
    discr = dir('*vhdr') ;
    for i=1:length(discr)
       EEG = pop_loadbv('.',discr(i).name) ; 
       EEG = pop_resample(EEG,256) ; 
       if i==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG) ; end
    end
    merged = pop_chanedit(merged,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
    merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61) - eegfiltfft(merged.data,merged.srate,84,86) - eegfiltfft(merged.data,merged.srate,0,1);  

    %figure,bar(sum(abs(diff(merged.data,1,2)),2))
    %merged = pop_interp(merged,[32],'spherical');
    pop_saveset(merged,'merged.set'); 
%}    
    
    
   
    merged = pop_loadset('merged.set');   
    for i=1:length(bades{sb})
       merged.data(bades{sb}(i),:) = rand(1,length(merged.data))/100;  
        
    end
   
    eps5 = pop_epoch(merged,{'S  5'},[-1,2]); 
    ep15 = pop_epoch(merged,{'S 15'},[-1,2]); 

    
end


