subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'};
for sb=1:length(subs)
    cd(['e:/badger_eeg/',subs{sb}]);
    ls 
    bcgweights= load('bcgweights'); 
    bcgweights = bcgweights.bcgweights; 
    if sb==6
        sets = {'grad_1*allstim*set','grad_3*allstim*set','grad_2*gamma*set','grad_4*gamma*set','grad_*movie*set','grad_*rest*set'};
    elseif sb==8
        sets = {'remove_gradient*allstims*01*set','remove_gradient*allstims*02*set','remove_gradient*gamma*01*set','remove_gradient*gamma*02*set','remove_gradient*movie*set','remove_gradient*rest*set'};
    else
        sets = {'grad*allstim*01*set','grad*allstim*02*set','grad*gamma*01*set','grad*gamma*02*set','grad*movie*set','grad*rest*set'};
    end
    
    for st=1:length(sets)
        if sb==5 && st==5
            setfreqs{5} = (setfreqs{1} + setfreqs{1} + setfreqs{1} + setfreqs{1}) / 4; 
        end
       %if st==1 ; merged = pop_loadset(sets(st).name); else merged = pop_mergeset(merged,pop_loadset(sets(st).name)); end
       if ~(sb==5 && st==5)
       eegname = dir(sets{st}); 
        eeg = pop_loadset(eegname.name); 
        bcgweights= load('bcgweights'); 
        bcgweights = bcgweights.bcgweights; 
        winv = pinv(bcgweights{1}*bcgweights{2}); 
        acts = (bcgweights{1}*bcgweights{2})*eeg.data; 
        
        trigtypes = {eeg.urevent.type}; 
        triglats = cell2mat({eeg.urevent.latency}); 
        r128s = find(strcmpi('R128',trigtypes)); 
        r128lats = triglats(r128s); 
        ntrs = length(r128lats); 
        freqs = 1:2:100; 
        allfreqs = zeros(50,5,length(acts)); 
        for f=1:length(freqs)
            freqsf = eegfiltfft(acts(1:5,:),eeg.srate,freqs(f)-1,freqs(f)+1);  
            %absf = mean(abs(freqsf),1); 
            absf = abs(freqsf); 
            allfreqs(f,:,:) = absf; 
        end
        allfreqs = allfreqs(:,:,r128lats(1):(length(r128lats)-1)*(eeg.srate*0.693)); 
        clear resfreqs; 
        for i=1:5   
            resfreqs(:,i,:) = imresize(squeeze(allfreqs(:,i,:)),[50,ntrs]); 
        end; 
        setfreqs{st} = resfreqs ;
       end
    end
    mkdir setfreqs ; cd setfreqs; 
    save('setfreqs','setfreqs'); 
    %{
    mergefilt = eegfiltfft(merged.data,merged.srate,1,128); 
    [weights,sphere] = runica(mergefilt(:,1:5:end),'maxsteps',128); 
    winv = pinv(weights*sphere); 
    figure,for i=1:64 ; subplot(5,13,i) ; topoplot(squeeze(winv(:,i)),merged.chanlocs);  title(i); end
    acts = winv*merged.data; 
    
    bcgweights = {weights,sphere};
    save('bcgweights','bcgweights'); 
    %}
end


