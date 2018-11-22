clear all ; close all 
cd c:/shared/fresh_eeg/

pupils = dir('anticipatory*vhdr') ; 
for p=1:length(pupils)
    EEG = pop_loadbv('.',pupils(p).name) ; 
    res = pop_resample(EEG,250) ; 
    if p==1 ; merged = res ; else merged = pop_mergeset(res,merged) ; end 
end
temp = merged.data(1:32,:) ; merged.data(1:32,:) = merged.data(33:64,:) ; merged.data(33:64,:) = temp ; 
mergefilt = merged ; mergefilt.data = eegfiltfft(merged.data,250,1,125) ; 
mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,250,59,61) ; 
mergefilt = pop_chanedit(mergefilt,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
mergesum = sum(diff(mergefilt.data,1,2).^2,2) ; bads = [13,48,54] ; 
mergefilt = pop_interp(mergefilt,bads,'spherical') ; 
epfilt = pop_epoch(mergefilt,{'S 11','S 12','S 13','S 14'},[-1,6]) ; 
mergeica = pop_runica(epfilt,'runica') ; 

applied = ica_applyweights(mergefilt,mergeica) ; 

tp(mergeica) ; 
[s,f] = spectopo(mergefilt.data,0,250,'plot','off') ; 

trigs = {'S 21','S 22','S 23','S 24'} ;
for i=1:length(trigs)
   epi = pop_epoch(applied,{trigs{i}},[-1,6]) ;  
   for j=1:64 ; 
      [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epi.icaact(j,:,:)),epi.pnts,[epi.xmin,epi.xmax],epi.srate,0,...
          'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0) ; 
   end
    
end
figure;for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(ersp(:,i,:,:),1)),[-3,3]) ; title(i);  end
tp(mergeica) ; 


trigs = {'S 21 (grating)','S 22 (low contrast)','S 23 (random)','S 24 (crosshair only)'} ;
comps = [3,4] ; 
for i=1:4 ; subplot(2,2,i) ; imagesc(times,freqs,squeeze(mean(ersp(i,comps,:,:),2)),[-5,5]) ;
    title(trigs{i}) ; vline(0,'r') ; vline(2,'k') ; vline(5,'k') ;axis xy ;
    if i==1 ; xlabel('time(s)') ; ylabel('freq(hz)') ; end
end

subplot(2,2,1) ; 
plot(squeeze(mean(mean(ersp(1,comps,27:35,:),2),3)),'r') ; hold on ; 
plot(squeeze(mean(mean(ersp(1,comps,5:12,:),2),3)),'b') ; hline(0,'k') ; ylim([-6,3]) ; vline([26],'r') ; vline([84,175],'k') ; 
set(gca,'XTick',1:25:length(times),'XTickLabel',round(times(1:25:end))) ; xlabel('time(s)') ; ylabel('power(db)') ; legend({'gamma','alpha/beta'}) ; title('grating') ; 
subplot(2,2,2) ; 
plot(squeeze(mean(mean(ersp(4,comps,27:35,:),2),3)),'r') ; hold on ; 
plot(squeeze(mean(mean(ersp(4,comps,5:12,:),2),3)),'b') ; hline(0,'k') ; ylim([-6,3]) ; vline([26],'r') ; vline([84,175],'k') ; 
set(gca,'XTick',1:25:length(times),'XTickLabel',round(times(1:25:end))) ; title('crosshair only') ; 



