

cd('\\snail\russell\honey_badger\alex\resamps') ; ls ; clear all ; close all ;
sounds=dir('r*set') ;
for i=1:max(size(sounds)) ; 
   %EEG = pop_loadbv('.',sounds(i).name) ; 
   %EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
   %EEG = pop_resample(EEG,256) ; 
   EEG = pop_loadset(sounds(i).name,'.') ; 
   if i==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end
    
end


zsums = zscore(sum(diff(merged.data,2).^2,2)) ; 
merged = pop_interp(merged,find(zsums>1),'spherical') ; 
ica = merged ; ica.data = eegfiltfft(ica.data,ica.srate,1,128) ;
ica = pop_runica(ica,'runica') ; 
trigs = {'S 10','S 11','S 12','S 14','S 15','S 13'} ;
clear ersp ;
for i=1:size(trigs,2) 
    ep = pop_epoch(ica,{trigs{i}},[-.85,2.85]) ;
    for j=1:64 ; 
        [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(j,:,1:20)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,'plotersp','off','plotitc','off',...
            'freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0) ; 
        
    end
end







