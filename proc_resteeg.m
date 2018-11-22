clear all ; close all 
cd c:/shared/fresh_eeg/
pupils = dir('rest*vhdr') ; 
for p=1:length(pupils)
    EEG = pop_loadbv('.',pupils(p).name) ; 
    res = pop_resample(EEG,250) ; 
    if p==1 ; merged = res ; else merged = pop_mergeset(res,merged) ; end 
end
temp = merged.data(1:32,:) ; merged.data(1:32,:) = merged.data(33:64,:) ; merged.data(33:64,:) = temp ; 
mergefilt = merged ; mergefilt.data = eegfiltfft(merged.data,250,1,35) ; 
mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,250,59,61) ; 
mergefilt = pop_chanedit(mergefilt,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
mergeica = pop_runica(mergefilt,'runica','maxsteps',128) ; 

acts = mergeica.icaweights*mergeica.icasphere*merged.data ;
[s,f] = spectopo(acts,0,mergeica.srate,'plot','off') ; 
tp(mergeica)

freqs = 1:50 ;
for fr=1:length(freqs)
  allfreqs(:,fr,:) = imfilter(abs(eegfiltfft(acts,mergeica.srate,freqs(fr)-1,freqs(fr)+1)),fspecial('gaussian',[1,300],50)) ;
end

for i=1:50
   corrmat(i,:,:) = corr(squeeze(allfreqs(:,i,:))') ;  
end

strings = {'1','2','3','4','5','6','7','8','9','10'} ; 
for i=1:50 ; subplot(5,10,i) ; imagesc(squeeze(corrmat(i,:,:)),[-1,1]) ; end









