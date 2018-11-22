subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie'} ;

for s=1%:length(subs) 
cd(['c:/shared/badger_eeg/',subs{s}]) ; 

eegs = dir('*gamma*Pulse*vhdr') ; 

for i=1:length(eegs)
    EEG = pop_loadbv('.',eegs(i).name) ; 
    labs = {EEG.urevent.type} ; 
    lats = cell2mat({EEG.urevent.latency}) ;
    r128s = find(strcmpi('R128',labs)) ; 
    st = lats(r128s(1)) ; en = lats(r128s(end)) ; 
    EEG = pop_select(EEG,'point',[st,en]) ; 
    if i==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end 
end
mergefilt = merged ; 
mergefilt.data = eegfiltfft(mergefilt.data,mergefilt.srate,1,128) ; 
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

trigs = {'S  1','S  2','S  3'} ; 
clear ersp ; 
for t=1:length(trigs)
    ep = pop_epoch(eegi,{trigs{t}},[-2,7]) ; 
    for c=1:64 ; 
        for trial=1:size(ep.icaact,3)
            [ersp(t,c,trial,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(c,:,trial)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                        'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'timesout',200) ; 
        end
    end
end

bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,200]) ; 
figure,
for i=1:64 ; subplot(5,13,i) ; imagesc(times,freqs,squeeze(mean(bersp(1,i,:,:,:),3)),[-6,6]) ; axis xy ; end 

save('bersp','bersp') ; eegi = pop_saveset(eegi,'eegi') ;
save('freqs','freqs') ; save('times','times') ; 
end

