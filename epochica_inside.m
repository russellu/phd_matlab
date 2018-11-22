clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie','tegan'} ; 
% inside
for s=1:length(subs)  ;
    cd(['c:/shared/badger_eeg/',subs{s}]) ;  
    prefix = '' ; 
    gammas=dir([prefix,'*gamma*Pulse*vhdr']) ;
    
    for g=1:length(gammas)
       if g==1 ; EEG = pop_loadbv('.',gammas(g).name) ; else EEG = pop_mergeset(pop_loadbv('.',gammas(g).name),EEG) ; end 
    end
    
    %trigs = {'S 11','S 12','S 13','S 14','S 21','S 22','S 23','S 24','S 31','S 32','S 33','S 34','S 41','S 42','S 43','S 44','S 51','S 52','S 53','S 54','S 61','S 62','S 63','S 64'} ; 
    trigs = {'S  1','S  2','S  3'} ; 
    epochmerged = pop_epoch(EEG,trigs,[-2,7]) ; 
    epochfilt = epochmerged ; 
    epochfilt.data = eegfiltfft(epochfilt.data,epochfilt.srate,45,70) ; 
    epochfiltica = pop_runica(epochfilt,'runica') ; 
    applied = ica_applyweights(EEG,epochfiltica) ; 
    neweeg = pop_epoch(applied,trigs,[-2,7]) ; 
    for i=1:64
          [ersplow(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(neweeg.icaact(i,:,:)),neweeg.pnts,[neweeg.xmin,neweeg.xmax],neweeg.srate,0,...
           'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',0,'timesout',100,'winsize',64) ;   
    end
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersplow(i,:,:))) ; title(i) ; end
    neweeg = pop_saveset(neweeg,'neweeg') ; 
end
