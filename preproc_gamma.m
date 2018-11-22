clear all ; close all
cd C:\shared\allres
subs = {'alex','alexandra3','audrey','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','leila','lisa','marc',...
    'marie','mathieu','maxime','menglu','mingham','patricia','po','russell','suhan2','sunachakan','tah','tegan2','vincent'} ; 
badchans2 = {[22,17],[28,12,16],[42,22],[53,26],[17,45,63],[16,12],[17,22],[28,22],[32,22,46],[22,33,41],[3,22,44],[12,46],[22],[32,28,12,2],[22,32,28],[17,35,34],[29,32,60],[22,55,42,12],[32,28],[22,53,48,42],[17],[17,2],[34,63,49],[22,28,2],[17,28],[32,17],[22,32,19,14,28],[32],[22,28,16,17],[32,25,64]} ; 
swaps = [1,1,0,0,0,0,0,1,1,1,1,0,0,1,1,1,0,0,1,1,0,0,0,1,0,0,1,0,0,1] ; 
trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'} ; 

for sb=23:length(subs)
    cd(['c:/shared/allres/',subs{sb}]) ; 
    ls 
    resamps = dir('resamp*set') ; 
    for res=1:length(resamps)
        eeg = pop_loadset(resamps(res).name) ;
        if swaps(sb) == 1 ; newdat = eeg.data(1:32,:) ; eeg.data(1:32,:) = eeg.data(33:64,:) ; eeg.data(33:64,:) = newdat ; end 
        zs(res,:) = zscore(sum(abs(diff(eeg.data,1,2)),2)) ; 
        if res==1 ; merged = eeg ; else merged = pop_mergeset(eeg,merged) ; end
    end
    merged = pop_interp(merged,badchans2{sb},'spherical') ; 
    merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61) ; 
    %pop_saveset(merged,'cleanfilt.set') ; 
    %oldeeg = pop_loadset('merged.set') ; 
    %weights = oldeeg.icaweights ; sphere = oldeeg.icasphere ; fullica{1} = weights ; fullica{2} = sphere ; 
    %save('highica','highica') ; 
    %{
    filtmerged = merged ; filtmerged.data = eegfiltfft(merged.data,merged.srate,40,100) ; 
    filtmerged = pop_epoch(filtmerged,trigs,[-1,2]) ; 
    filtmerged = pop_runica(filtmerged,'runica','maxsteps',128) ;    
    newmerged = merged ; newmerged.data = filtmerged.icaweights*filtmerged.icasphere*merged.data ; 
    epica = pop_epoch(newmerged,{'S 11','S 13','S 15'},[-1,2.5]) ; clear ersp ;
    for i=1:64 ; disp(i)
        [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(i,:,:)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                            'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',0,'verbose','off','timesout',200) ; 
    end
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(i,:,:)),[-5,5]) ; end ; suptitle(subs{sb}) ; 
    %}
    merged = pop_resample(merged,256) ; pop_saveset(merged,'cleanfilt.set') ; 
    %highepochica{1} = filtmerged.icaweights ; highepochica{2} = filtmerged.icasphere ; save('highepochica','highepochica') ; 
end
