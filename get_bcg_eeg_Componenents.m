clear all ; close all; 
subs = {'alex','dina','genevieve','jeremie','russell','tegan','valerie'}; 
sets = {'retino_allstim*01*set','retino_allstim*02*set','retino_gamma*01*set','retino_gamma*02*set','retino_movie*set','retino_rest*set'}; 

freqs = {[1,3],[5,15],[16,25],[35,60]}; 

for sb=1:length(subs)
        
    % get the bcg (and electrode bcg) 
    for st=1:length(sets)
       cd(['E:\badger_eeg\',subs{sb}]);
       set_st = dir(sets{st}); 
       eeg = pop_loadset(set_st.name); 
       
       if st==1 ; merged = pop_loadset(set_st.name) ; else merged = pop_mergeset(merged,pop_loadset(set_st.name)); end
    end
    
    mergefilt = eegfiltfft(merged.data,merged.srate,1,128); 
    [weights,sphere] = runica(mergefilt(:,1:2:end),'maxsteps',128); 
    winv = pinv(weights*sphere); 
    acts = weights*sphere*mergefilt; 
    figure,
    for i=1:10
        subplot(4,6,i*2-1) ; topoplot(winv(:,i),merged.chanlocs);
        subplot(4,6,i*2); plot(acts(i,100000:105000)); 
    end
    
    bcgica{1} = weights; bcgica{2} = sphere ; save('bcgica','bcgica'); 
end






