clear all ; close all  ;
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','valerie','tegan'} ;
bands = {[45,60],[50,65],[50,70],[55,70],[50,70],[50,70],[30,70],[50,70],[45,70]} ; 

for s=9:length(subs) 
cd(['c:/shared/badger_eeg/',subs{s}]) ; 
eegs = dir('*Pulse*vhdr') ; 
clear alleegs ; 
    for i=1:length(eegs)
        EEG = pop_loadbv('.',eegs(i).name) ; 
        labs = {EEG.urevent.type} ; 
        lats = cell2mat({EEG.urevent.latency}) ;
        r128s = find(strcmpi('R128',labs)) ; 
        st = lats(r128s(1)) ; en = lats(r128s(end)) ; 
        raweegs{i} = EEG ; % apply the weights to this dataset later
        EEG = pop_select(EEG,'nopoint',[1,st ; en,size(EEG.data,2)]) ; 
        % figure,plot(std(EEG.data([1:31,33:64],:),0,1))
        stdchans = std(EEG.data([1:31,33:64],:),0,1) ; 
        zstdchans = zscore(stdchans) ; bads = find(zstdchans>3) ; 
        badwindow = EEG.srate*1.5 ; % +/- 1.5 second
        clear badinds
        for b=1:length(bads)
            badinds(b,:) = bads(b)-badwindow:bads(b)+badwindow ;         
        end
        badinds(badinds==0) = 1 ; 
        uniquebads = unique(badinds) ;
        uniquebads(uniquebads<0 | uniquebads > size(EEG.data,2)) = []  ; 
        zclust = zeros(1,size(EEG.data,2)) ; zclust(uniquebads) = 1 ; 
        bw = bwconncomp(zclust) ; 
        plist = bw.PixelIdxList ; 
        clear borders
        for bd=1:length(plist)
            borders(bd,:) = [plist{bd}(1),plist{bd}(length(plist{bd}))] ; 
        end
        EEG = pop_select(EEG,'nopoint',borders) ; 
        stdchans2 = std(EEG.data([1:31,33:64],:),0,1) ; 
        %figure,plot(stdchans) ; hold on ; plot(stdchans2,'r') ; 
        %figure, plot(mat2gray(stdchans)) ; hold on ; plot(zclust,'r') ; 
        if i==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end 
        alleegs{i} = EEG ; 
    end   
    mergefilt = merged ; 
    mergefilt.data = eegfiltfft(mergefilt.data,mergefilt.srate,bands{s}(1),bands{s}(2)) ; % bands{s}(1),bands{s}(2)
    mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,mergefilt.srate,59.5,60.5) ; 
    ica = pop_runica(mergefilt,'runica') ; 
    for i=1:length(alleegs)
        eegi = raweegs{i} ; 
        eegi.icaact = icaact(eegi.data,ica.icaweights*ica.icasphere,0) ; 
        eegi.icachansind = ica.icachansind ; 
        eegi.icasphere = ica.icasphere ; 
        eegi.icasplinefile = ica.icasplinefile ; 
        eegi.icaweights = ica.icaweights ; 
        eegi.icawinv = ica.icawinv ; 
        eegi.data = eegi.data - eegfiltfft(eegi.data,mergefilt.srate,59.5,60.5) ; 
        eegi.icaact = eegi.icaact - eegfiltfft(eegi.icaact,mergefilt.srate,59.5,60.5) ; 
        eegi = pop_saveset(eegi,['highfreq_preproc_',strrep(eegs(i).name,'.vhdr','')]) ; 
    end
end
