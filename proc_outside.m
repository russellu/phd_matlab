clear all  ; close all
subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie'} ;

for sub=1:length(subs) ; 
cd(['c:/shared/badger_eeg/',subs{sub},'/outside']) ; ls  ; 
sounds=dir('*outside*vhdr') ;
for i=1:max(size(sounds)) ; 
   EEG = pop_loadbv('.',sounds(i).name) ; 
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
   EEG = pop_resample(EEG,250) ;
   eegs{i} = EEG ; 
   if i==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end   
end
mergefilt = merged ; mergefilt.data = eegfiltfft(mergefilt.data,mergefilt.srate,1,128) ; 
mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,mergefilt.srate,59.5,60.5) ; 
ica = pop_runica(mergefilt,'runica') ; 
eegi = merged ; 
eegi.icaact = icaact(eegi.data,ica.icaweights*ica.icasphere,0) ; 
eegi.icachansind = ica.icachansind ; 
eegi.icasphere = ica.icasphere ; 
eegi.icasplinefile = ica.icasplinefile ; 
eegi.icaweights = ica.icaweights ; 
eegi.icawinv = ica.icawinv ; 
eegi.data = eegi.data - eegfiltfft(eegi.data,mergefilt.srate,59.5,60.5) ; 
eegi.icaact = eegi.icaact - eegfiltfft(eegi.icaact,mergefilt.srate,59.5,60.5) ; 
comps = 1:64 ;
trigs = {'S  1','S  2','S  3','S  4'} ;
clear ersp 
for trig=1:length(trigs) ; 
    ep = pop_epoch(eegi,{trigs{trig}},[-1.8,5.5]) ; 
    for c=1:length(comps) ; 
        for tr=1:size(ep.icaact,3)
            [ersp(trig,c,tr,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(comps(c),:,tr)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'timesout',200) ; 
        end
    end
end
bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,200]) ; 
for t=2 ; figure
    for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(bersp(t,i,:,:,:),3)),[-6,6]) ; title(i) ; end
end
save('bersp','bersp') ; eegi = pop_saveset(eegi,'eegi') ;
save('freqs','freqs') ; save('times','times') ; 

end