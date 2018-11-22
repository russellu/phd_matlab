clear all ; close all ; 
cd c:/shared/badger_eeg/alex ; ls 
dens = dir('bcg_grad*gamma*set') ; 
bvs = dir('*gamma*Pulse*set') ; 

for i=1:length(dens) ;
   EEG = pop_loadset(dens(i).name) ; 
    if i==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end ; 
   
     EEG2 = pop_loadset(bvs(i).name) ; 
    if i==1 ; merged2 = EEG2 ; else merged2 = pop_mergeset(EEG2,merged2) ; end ; 
    
   %{
   etypes = {EEG.urevent.type} ; 
   st = find(strcmp('S 98',etypes)) ; 
   lats = cell2mat({EEG.urevent.latency}) ;
   stlat = lats(st) ; 
   endlat = stlat + EEG.srate * (8.5*60) ;  
   if endlat > size(EEG.data,2) ; endlat = size(EEG.data,2) ; end
   EEG = pop_select(EEG,'point',[stlat,endlat]) ;   
  %}
end
for i=1:2 
    if i==1 ; FILT = merged ; else FILT = merged2 ; end
    FILT.data = eegfiltfft(FILT.data,FILT.srate,1,60) ; 
    EP = pop_epoch(FILT,{'S  1','S  2','S  3'},[-2,8]) ; 
    ica = pop_runica(EP,'runica') ; 
    ica = pop_chanedit(ica,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 

    if i==1; newmerged = ica_applyweights(merged,ica) ; ep = pop_epoch(newmerged,{'S  1','S  2','S  3'},[-2,7]) ; 
    else newmerged = ica_applyweights(merged2,ica) ; ep = pop_epoch(newmerged,{'S  1','S  2','S  3'},[-2,7]) ; end
   
    clear ersp ; 
    for c=1:64 ;
       [ersp(c,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(c,:,:)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
           'plotersp','off','plotitc','off','timesout',100,'freqs',[1,120],'nfreqs',60,'winsize',64) ;      

    end
    figure;for j=1:64 ; subplot(5,13,j) ; imagesc(squeeze(ersp(j,:,:)),[-5,5]) ; end

end




