cd('C:\shared\badger\jeremie') ; ls  ; clear all  ; close all
sounds=dir('*retino_*vhdr') ;
for i=1:max(size(sounds)) ; 
   EEG = pop_loadbv('.',sounds(i).name) ; 
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
   EEG = denoise_grad3(EEG) ; 
   EEG = pop_resample(EEG,256) ;
   EEG = denoise_bcg2(EEG) ; 
   eegs{i} = EEG ; 
   if i==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end   
end

for i=1:5 
      if i==1 ; merged = eegs{i} ; else merged = pop_mergeset(eegs{i},merged) ; end  
end
mergefilt = merged ; mergefilt.data = eegfiltfft(merged.data,merged.srate,10,80) ; 
ica = pop_runica(mergefilt,'runica') ;
eeg2 = merged ; 
eeg2.icaact = icaact(eeg2.data,ica.icaweights*ica.icasphere,0) ; 
eeg2.icachansind = ica.icachansind ; 
eeg2.icasphere = ica.icasphere ; 
eeg2.icasplinefile = ica.icasplinefile ; 
eeg2.icaweights = ica.icaweights ; 
eeg2.icawinv = ica.icawinv ;   
ep = pop_epoch(eeg2,{'S  1','S  2'},[-2,7]) ; 
comps=1:64 ; clear ersp ; 
for c=1:length(comps) ; 
    for t=1:size(ep.icaact,3)
        [ersp(c,t,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(comps(c),:,t)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                'plotersp','off','plotitc','off','freqs',[1,128],'nfreqs',60,'winsize',64,'baseline',NaN,'timesout',200) ;     
    end
end

bersp = ersp - repmat(mean(ersp(:,:,:,times<0),4),[1,1,1,200]) ; 
mbersp = squeeze(mean(bersp(:,1:64,:,:),2)) ;

for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mbersp(i,:,:)),[-4,4]) ; title(i) ; end


