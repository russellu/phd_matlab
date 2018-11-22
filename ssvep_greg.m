cd C:\shared\greg_ssvep ; ls 

ssveps = dir('*vhdr') ; 
for s=1:length(ssveps)
g1 = pop_loadbv('.',ssveps(s).name) ; 
g1 = pop_resample(g1,256) ; 
if s>1 ; merged = pop_mergeset(g1,merged) ; else merged = g1 ; end
end

mergefilt = merged ;
mergefilt.data = eegfiltfft(merged.data,merged.srate,15,17) ; 
ep = pop_epoch(mergefilt,{'S 10'},[0,60]) ; 
ica = pop_runica(ep,'runica') ; 
ica = pop_chanedit(ica,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 

newmerged = ica_applyweights(merged,ica) ; 
newmerged.data = eegfiltfft(newmerged.data,newmerged.srate,2,128) ; 
ep = pop_epoch(newmerged,{'S 10'},[0,60]) ; 
[s,f] = spectopo(squeeze(mean(ep.icaact(:,:,:),3)),0,256) ;



clear ersp ; 
for c=1:64 ;
   [ersp(c,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(c,:,:)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
       'plotersp','off','plotitc','off','timesout',1000,'freqs',[1,120],'nfreqs',120,'winsize',256) ;      
end

plot(squeeze(mean(mean(ersp(:,freqs>=6 & freqs<=8,:),2))))



