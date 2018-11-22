cd c:/shared/greg ; ls 

g1 = pop_loadbv('.','greg_02.vhdr') ; 
g1 = pop_resample(g1,256) ; 
filtg1 = g1 ; 
filtg1.data = eegfiltfft(g1.data,g1.srate,60,80) ; 

ep = pop_epoch(filtg1,{'S 11','S 12','S 13','S 14'},[-1,5]) ; 
icag1 = pop_runica(ep,'runica') ; 
icag1 = pop_chanedit(icag1,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 


g1 = ica_applyweights(g1,icag1) ; 
stims = {'S 11','S 12','S 13','S 14'} ; 
clear ersp ; 
for i=1:4
    ep = pop_epoch(g1,{stims{i}},[-1,5]) ; 
    for c=1:64 ;
       [ersp(i,c,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(c,:,:)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
           'plotersp','off','plotitc','off','timesout',200,'freqs',[1,120],'nfreqs',60,'winsize',64) ;      
    end
    
end

for i=1:4 ; figure ; for j=1:64 ; subplot(5,13,j) ; imagesc(squeeze(ersp(i,j,:,:)),[-4,4]) ; end ; end