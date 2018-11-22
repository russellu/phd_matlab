clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie'} ; 
eegcomps = {[30,26],[32,23],[21,28],[47,46],[19,37],[47,38]} ;

for sub=1:length(subs)
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
    allstims=dir([prefix,'preproc*allstim*set']) ;
    setNames = {allstims(1).name,allstims(2).name} ;
    
    etrigs = {'S 11','S 12','S 13','S 14','S 21','S 22','S 23','S 24','S 31','S 32','S 33','S 34','S 41','S 42','S 43','S 44',...
              'S 51','S 52','S 53','S 54','S 61','S 62','S 63','S 64','S 71','S 72','S 73','S 74','S 81','S 82','S 83','S 84'} ; 
    
    clear trialdata trialTypeLats allEventIndices allbold allersp ; 
    for s=1:length(setNames)
        EEG = pop_loadset(setNames{s}) ; 
        %EEG.icaact = eegfiltfft(EEG.icaact,EEG.srate,1,120) ; 
        r128s = find(strcmp('R128',{EEG.urevent.type})) ; 
        lats = {EEG.urevent.latency} ; 
        trlats = cell2mat(lats(r128s)) ; 
        offsetSec = trlats(1)./EEG.srate ; 
        ep = pop_epoch(EEG,etrigs,[-1,6]) ; 
        ep2 = pop_epoch(EEG,etrigs,[-8,16]) ; 

        types = {ep.epoch.eventtype} ; 
        
        allEventIndices = {ep.epoch.eventurevent} ; 
        eventTypes = {ep.epoch.eventtype} ; 
        trialTypeLats = zeros(1,32) ; 
        trialdata = zeros(32,size(ep2.icaact,1),size(ep2.icaact,2)) ;  
        for i=1:length(allEventIndices)
            eventIndicesI = allEventIndices{i} ; 
            for e=1:length(etrigs)
                ind = find(strcmp(etrigs{e},eventTypes{i})) ; 
                if ~isempty(ind) ; break ; end
            end
            lat_i = lats{eventIndicesI{ind}} ;      
            trialdata(sum(trialTypeLats~=0)+1,:,:) = ep2.icaact(:,:,i) ; 
            trialTypeLats(sum(trialTypeLats~=0)+1) = lat_i ; 
        end
        
        trialTypeLats = trialTypeLats - trlats(1) ; 
        clear ersp 
        for comp=1:length(eegcomps{sub}) ; 
            for i=1:size(trialdata,1)
                    [ersp(comp,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(trialdata(i,eegcomps{sub}(comp),:)),ep2.pnts,[ep2.xmin,ep2.xmax],ep2.srate,0,...
                                                                            'plotersp','off','plotitc','off','baseline',NaN,'timesout',100,'freqs',[1,100],'nfreqs',50,'winsize',64) ; 
            end
        end
               
        allersp(:,(s-1)*32+1:(s-1)*32+32,:,:) = ersp ; 
        
        smix = segmix{s} ; 
        smix = eegfiltfft(smix',1/0.693,0.01,1.5) ; smix = smix' ; 
        trs = round((reshape(trialTypeLats,[1,numel(trialTypeLats)])./EEG.srate)./0.693) ; 
        stimtrs = round(trialTypeLats./(EEG.srate*0.693)) ; 
        design = zeros(1,size(smix,1)) ; stimT = round(10./0.693) ; 
        for i=1:length(trs) ; design(trs(i):trs(i)+stimT) = 1 ; end 
        hrf = spm_hrf(0.693) ; 
        conved = conv(design,hrf,'full') ; 
        conved = conved(1:size(smix,1)) ; 
        
        corrs = corr(smix(50:end-50,:),conved(50:end-50)') ; 
        [sv,si] = sort(corrs,'descend') ; 
        topts = mean(smix(:,si(1)),2) ; 
                
        baseT = round(1./0.693) ; taskT = round(15./0.693) ; 
        clear boldepochs
        for j=1:size(stimtrs,2)
            boldepochs(j,:) = topts(stimtrs(j)-baseT:stimtrs(j)+taskT) - squeeze(mean(topts(stimtrs(j)-baseT))) ;            
        end
        allbold((s-1)*32+1:(s-1)*32+32,:) = boldepochs ; 
    end     
    
    mersp = squeeze(mean(allersp,1)) ; 
    for c=1:2
        for i=1:50 ;
           for j=1:100
               fcorrs(c,i,j) = corr2(mean(allbold(:,16:20),2),squeeze(allersp(c,:,i,j))') ;      
           end       
        end
    end
    allfcorrs(sub,:,:,:) = fcorrs ; 
end










