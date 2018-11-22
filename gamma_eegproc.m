clear all ; close all
cd C:\shared\allres
subs = {'alex','alexandra3','audrey','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','leila','lisa','marc',...
    'marie','mathieu','maxime','menglu','mingham','olga','patricia','po','russell','suhan2','sunachakan','tah','tegan2','vincent'} ; 
badchans = {[49,54],[35,44,48],[22,42],[42],[7,17,45],[12,16,21],[17],[54,60],[54,64],[1,54],[12,35,54],[42],[22],[60,64],[54,64],[8,49],[45],[22],[13,43],...
    [16,21,54],[17],[42,45],[10,11,35,50],[34,63],[54],[17],[17,32],[51,54,60,64],[32],[8,22,42],[64]} ;
swaps = [1,1,0,0,0,0,0,1,1,1,1,0,0,1,1,1,0,0,1,1,0,0,1,0,1,0,0,1,0,0,1] ; 

for sb=1:length(subs)
    cd(['c:/shared/allres/',subs{sb}]) ; 
    ls 
    resamps = dir('res*set') ; 
    for rs=1:length(resamps)
       eeg = pop_loadset(resamps(rs).name) ; 
       if rs==1 ; merged = eeg ; else merged = pop_mergeset(eeg,merged) ; end     
    end
    merged = pop_interp(merged,badchans{sb},'spherical') ; 
    if swaps(sb) == 1 ; tempdat = merged.data(1:32,:) ; merged.data(1:32,:) = merged.data(33:64,:) ; merged.data(33:64,:) = tempdat ; end 
    merged = pop_resample(merged,256) ; 
    pop_saveset(merged,'merged.set') ; 
    trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'} ; 
    mergefilt = merged ; mergefilt.data = eegfiltfft(merged.data,merged.srate,50,100) ;%+ eegfiltfft(merged.data,merged.srate,10,15)/5 ; 
    eps = pop_epoch(mergefilt,trigs,[-1,3]) ; 
    epica = pop_runica(eps,'runica','maxsteps',128) ; 
    weights = epica.icaweights ; sphere = epica.icasphere ; 
    winv = pinv(weights*sphere) ; 
    figure,for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),merged.chanlocs) ; end ; suptitle(subs{sb}) ; 
    trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'} ; 
    newmerged = merged ; newmerged.data = weights*sphere*merged.data ; 
    allep = pop_epoch(newmerged,{trigs{1},trigs{5}},[-1,3]) ; 
    for i=1:size(allep.data,1)
        [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',200) ; 
    end
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(i,:,:)),[-4,4]) ; title(i) ; end ; suptitle(subs{sb}) ; 
    hhwcomps{1} = weights ; hhwcomps{2} = sphere ; save('hhwcomps','hhwcomps') ; 

end



