clear all ; close all ; 
cd C:\shared\retino_eeg
discr = dir('russ_ret*vhdr') ; 
for i=1:length(discr)
   EEG = pop_loadbv('.',discr(i).name) ; 
   EEG = pop_resample(EEG,256) ; 
   if i==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG) ; end
end
merged = pop_chanedit(merged,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
temp = merged.data(1:32,:) ; merged.data(1:32,:) = merged.data(33:64,:) ; merged.data(33:64,:) = temp ; 
mergefilt = merged ; mergefilt.data = eegfiltfft(merged.data,merged.srate,40,90) ; 
ica = pop_runica(mergefilt,'runica') ; ica.icaact = ica.icaweights*ica.icasphere*merged.data ; 
ica.data = ica.icaact ; 
allep = pop_epoch(ica,{'S  1'},[0,120]) ;
clear allersp ; 
for i=1:64 ;
        [allersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',512,'baseline',0,'verbose','off') ; 
end
bersp = allersp - repmat(mean(allersp,3),[1,1,200]) ; 
for i=1:64 ; subplot(5,13,i) ; imagesc(medfilt2(squeeze(bersp(i,:,:)),[5,5]),[-2,2]) ; title(i) ; axis xy ; end




