clear all ; close all  ;
cd c:/shared/resmerged ; 
subs=dir('*') ; subs(1:2) = [] ; 


for s=20:length(subs) ; 
    cd(['c:/shared/resmerged/',subs(s).name]) ;  
    ls 
    merged = pop_loadset('merged.set') ; 
    mergefilt = eegfiltfft(merged.data,merged.srate,6,90) ; merged.data = mergefilt ; 
    ep = pop_epoch(merged,{'S 11','S 12','S 13','S 14','S 15','S 16'},[-1,3]) ; 
    mergeica = pop_runica(ep,'runica') ;
    filtica{1} = mergeica.icaweights ; filtica{2} = mergeica.icasphere ; 
    save('filtica','filtica') ; 
end