clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie','tegan','sukhman'} ; 
for sb=1:length(subs) 
    sub = subs{sb} ; 
    cd(['c:/shared/badger_eeg2/',sub]) 
    freqs = 1:2:60 ; 
    grads = dir('remove_grad*set') ;
    for g=1:length(grads)
        EEG = pop_loadset(grads(g).name) ;    
        [zdat,e4,ss] = fbcg(EEG) ; 
        allsets{g} = e4 ; 
        if g==1 ; merged = e4 ; else merged = pop_mergeset(e4,merged) ; end
    end
    merged = pop_chanedit(merged,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
    stds = std(merged.data,0,1) ; smoothstds = smooth(stds.^2,100) ; 
    [sv,si] = sort(stds,'descend') ;
    newdat = merged.data ; newdat(:,si(1:20000)) = [] ; 
    ica = pop_runica(merged,'runica') ; 
    for g=1:length(grads)
        name = strrep(grads(g).name,'remove_gradient2','bcgica') ;
        eeg = allsets{g} ; 
        eeg.icasphere = ica.icasphere ; eeg.icaweights = ica.icaweights ; 
        pop_saveset(eeg,name) ; 
    end
end








