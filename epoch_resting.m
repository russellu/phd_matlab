clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie'} ; 
eegcomps = {[30,26],[32,23],[21,28],[47,46],[19,37],[47,38]} ;

for sub=1%:length(subs)
    %%%% BOLD processing:
    cd(['c:/shared/badger_mri/',subs{sub},'/nii/melodic']) ; 
    mix = load('melodic_mix') ; 
    weights = load_untouch_nii('melodic_IC.nii.gz') ; 
    allweights{sub} = weights.img ; 
    meansub = load_untouch_nii('mean.nii.gz') ; 
    allmeans{sub} = meansub.img ; 
    mixinds = 1:735:735*5 ; 
    for i=1:length(mixinds) ; boldinds{i} = mixinds(i):mixinds(i)+734 ; end
    boldinds{6} = 735*5+1:735*5+1+449 ; 
    for i=1:length(boldinds) ; segmix{i} = mix(boldinds{i},:) ; end

    %%%% EEG processing:
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    prefix = 'highfreq_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ; %gammas2=dir([prefix,'preproc*movie*set']) ;
    setNames = {gammas1(2).name} ;
    
    etrigs = {'R128'} ; 
    
    clear trialdata trialTypeLats allEventIndices allbold allersp ; 
    for s=1%:length(setNames)
        EEG = pop_loadset(setNames{s}) ; 
        EEG.icaact = eegfiltfft(EEG.icaact,EEG.srate,1,120) ; 
        r128s = find(strcmp('R128',{EEG.urevent.type})) ; 
        lats = {EEG.urevent.latency} ; 
        types = {EEG.urevent.type} ; 
        trlats = cell2mat(lats(r128s)) ; 
        offsetSec = trlats(1)./EEG.srate ; 
        ep = pop_epoch(EEG,{'R128'},[-0.1,0.1]) ; % get all the triggers (to find times)
        
        evts = {ep.epoch.eventtype} ; 
        allEventIndices = {ep.epoch.eventurevent} ; 
        clear latinds 
        for i=1:length(evts)
            evti = evts(i) ; 
            currentind = find(strcmp('R128',evti{:})) ; 
            latinds(i) = allEventIndices{i}(currentind) ; 
        end
        latinds = cell2mat(latinds) ;
        latinds([1:15,length(latinds)-15:length(latinds)]) = [] ; % remove start and end of scan effects
        ep2 = pop_epoch(EEG,{},[-7,7],'eventindices',latinds) ; % get epochs only within a range 
        
        clear ersp 
        for comp=1:length(eegcomps{sub}) ; 
            for i=1:size(ep2.icaact,3)
                [ersp(comp,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep2.icaact(eegcomps{sub}(comp),:,i)),ep2.pnts,[ep2.xmin,ep2.xmax],ep2.srate,0,...
                                                                        'plotersp','off','plotitc','off','baseline',NaN,'timesout',100,'freqs',[1,100],'nfreqs',50,'winsize',64) ; 
            end
        end
        
        fmix = segmix{4} ; 
        fmix = eegfiltfft(fmix',1/0.693,0.01,1.5) ; fmix = fmix' ; 
        fmix = fmix(16:end-16,:) ; 
        mersp = squeeze(mean(ersp,1)) ; 
        clear kernel
        for c=1:85
            for i=1:50 
                for j=1:100
                    kernel(c,i,j) = corr2(squeeze(mersp(:,i,j)),fmix(:,c)) ; 
                end
            end
        end
    end        
end








