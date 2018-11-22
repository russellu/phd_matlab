clear all ; close all 
cd c:/shared/fresh_eeg/
trigs = {'S 12','S 13','S 14','S 15','S 32','S 33','S 34','S 35'} ;
pupils = dir('stimcomps*vhdr') ; 
for p=1:length(pupils)
    EEG = pop_loadbv('.',pupils(p).name) ; 
    res = pop_resample(EEG,250) ; 
    if p==1 ; merged = res ; else merged = pop_mergeset(res,merged) ; end 
end
temp = merged.data(1:32,:) ; merged.data(1:32,:) = merged.data(33:64,:) ; merged.data(33:64,:) = temp ; 
mergefilt = merged ; mergefilt.data = eegfiltfft(merged.data,250,1,125) ; 
mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,250,59,61) ; 
mergefilt = pop_chanedit(mergefilt,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
%mergesum = sum(diff(mergefilt.data,1,2).^2,2) ; bads = [13,48,54] ; 
%mergefilt = pop_interp(mergefilt,bads,'spherical') ; 
epfilt = pop_epoch(mergefilt,trigs,[-1.5,1.5]) ; 
mergeica = pop_runica(epfilt,'runica') ; 

applied = ica_applyweights(mergefilt,mergeica) ; 

tp(mergeica) ; 
[s,f] = spectopo(mergefilt.data,0,250,'plot','off') ; 

clear s s2 ; 
for i=1:length(trigs)
   epi = pop_epoch(applied,{trigs{i}},[-.5,1.5]) ;  
   for j=1:64 ;
      [s(i,j,:,:),f] = spectopo(squeeze(epi.icaact(j,epi.times<1000 & epi.times>0,:))',0,epi.srate,'plot','off') ;
      [s2(i,j,:,:),f2] = spectopo(squeeze(epi.icaact(j,epi.times<0,:))',0,epi.srate,'plot','off'); 
   end
   
   %{
   for j=1:64 ; 
      [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epi.icaact(j,:,:)),epi.pnts,[epi.xmin,epi.xmax],epi.srate,0,...
          'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0) ; 
   end
    %}
end
clear meanstim
for freq=1:100 ; 
meanstim(freq,:,:,:) = mean(s(:,:,:,f>freq-1 & f<freq+1),4) - mean(s2(:,:,:,f2>freq-1 & f2<freq+1),4) ; 
end

bars = mean(meanstim,4)./std(meanstim,0,4) ; 





figure,for i=1:8 ; subplot(2,4,i) ; imagesc(squeeze(mean(ersp(i,comps,:,:),2)),[-6,6]) ; end

figure;for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(ersp(:,i,:,:),1)),[-3,3]) ; title(i);  end
tp(mergeica) ; 
comps = [8] ;





[s,f] = spectopo(mergeica.icaact,0,mergeica.srate,'plot','off') ; 



