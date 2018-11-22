clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 
eegcomps = {[30,26],[32,23],[21,28],[47,46],[19,37],[47,38],[44,46]} ;


for sub=1:length(subs)
    
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    prefix = 'highfreq_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ; gammas2=dir([prefix,'preproc*gamma*set']) ;
    setNames = {gammas1(1).name,gammas2(2).name} ;
   
    etrigs = {'S  1','S  2','S  3'} ; 
  
    chans = zeros(1,64) ; chans(eegcomps{sub}) = 1 ; bads = find(chans==0) ; 
    
    for s=1:length(setNames)
        EEG = pop_loadset(setNames{s}) ; 
        EEG = pop_subcomp(EEG,bads,0) ; 
        if s==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end
    end
    
    trials = pop_epoch(merged,{'S  1'},[-5,8]) ; 
    elabs = {EEG.chanlocs.labels} ;
     for comp=1:size(trials.data,1) ; 
         [ersp(comp,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(trials.data(comp,:,:)),trials.pnts,[trials.xmin,trials.xmax],trials.srate,0,...
                                                                        'plotersp','off','plotitc','off','baseline',NaN,'timesout',100,'freqs',[1,100],'nfreqs',50,'winsize',64) ; 
     end
     meanfreqs(sub,:,:) = squeeze(mean(ersp(:,:,times>0 & times<5),3)) ; 
     
     for comp=1:size(trials.data,1) ; 
         [ersp2(comp,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(trials.data(comp,:,:)),trials.pnts,[trials.xmin,trials.xmax],trials.srate,0,...
                                                                        'plotersp','off','plotitc','off','baseline',0,'timesout',100,'freqs',[1,100],'nfreqs',50,'winsize',64) ; 
     end
     meanfreqsbase(sub,:,:) = squeeze(mean(ersp2(:,:,times>0 & times<5),3)) ; 

end





for s=1:length(subs)
   cd(['c:/shared/eeg_comparisons/',subs{s}]) ;  
   freqpower = squeeze(meanfreqs(s,:,:)) ; 
   save('freqpower','freqpower') ; 
   save('freqs','freqs') ; 
end





