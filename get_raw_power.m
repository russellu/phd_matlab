clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 

comps = {[5,11,14],[5,9,12],[12,21],[7,9,27],[4,5],[4,9,11,14],[3,7,12],[5,14,16],[1,14],[8,12,30],[3,6,15,19,27],[2,6,11,20],...
         [3,21],[2,4,7],[5,7,13],[13,19],[5,7,11],[3,6,12],[5,9,12],[4,7],[4,7,13],[1,3,5],[8,12,17],[2,4,8]};


for sb=1:length(subs)
    disp(sb); 
    cd(['E:\clean_allres\',subs{sb}]) ; ls 
    
    eeg = pop_loadset('interp_flip_eeg.set'); 
    cleanweights = load('cleanweights'); cleanweights = cleanweights.cleanweights; 
    weights = cleanweights{1}; sphere = cleanweights{2}; 
    winv = pinv(weights*sphere); 
    neweeg = eeg;
    acts = weights*sphere*eeg.data;     
    sortcomps = load('sortcomps') ; 
    sortcomps = sortcomps.sortcomps; 
    invacts = winv(:,sortcomps(1:5))*acts(sortcomps(1:5),:); 
    neweeg.data = invacts; 
    
    freqs = 1:2:120;
    clear freqtopos; 
    invacts = invacts - eegfiltfft(invacts,eeg.srate,59,61) - eegfiltfft(invacts,eeg.srate,84,86); 
    for f=1:length(freqs)    
        gammafilt = eegfiltfft(invacts,eeg.srate,freqs(f)-1,freqs(f)+1) ;
        gammaeeg = eeg; gammaeeg.data = gammafilt; 
        pop_saveset(gammaeeg,['gen_hz_',num2str(f*2),'.set']);    
    end
    
    gammafilt = eegfiltfft(invacts,eeg.srate,8,25) ;
    gammaeeg = eeg; gammaeeg.data = gammafilt; 
    pop_saveset(gammaeeg,['gen_alpha_hz.set']);   
    gammafilt = eegfiltfft(invacts,eeg.srate,40,90) ;
    gammaeeg = eeg; gammaeeg.data = gammafilt; 
    pop_saveset(gammaeeg,['gen_gamma_hz.set']);   
       
end


