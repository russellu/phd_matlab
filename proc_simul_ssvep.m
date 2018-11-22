cd C:\shared\Greg_2016-07-28\greg_simultaneous
ssvep = dir('*Pulse*vhdr') ; 
for s=1:length(ssvep)
    EEG = pop_loadbv('.',ssvep(s).name) ; 
    if s==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end
end

filt = merged ;
filt.data = eegfiltfft(filt.data,filt.srate,1,128) ; 
ep = pop_epoch(filt,{'S 11'},[0,240]) ; 
ep = pop_chanedit(ep,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
ica = pop_runica(ep,'runica') ; 

newmerged = ica_applyweights(merged,ica) ; 
newepoch = pop_epoch(newmerged,{'S 11'},[0,240]) ;

m1 = squeeze(mean(newepoch.icaact(:,:,1:2:end),3)) ; 

[s,f] = spectopo(newepoch.icaact,0,250,'plot','off') ; 



















