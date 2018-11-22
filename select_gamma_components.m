clear all ; close all ; 
cd C:\shared\allres
subs = {'alex','alexandra3','audrey','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','leila','lisa','marc',...
    'marie','mathieu','maxime','menglu','mingham','olga','patricia','po','russell','suhan2','sunachakan','tah','tegan2','vincent'} ; 
trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'} ; 

comps = {[4,8,12,14,21,22,40],[3,6,7,9,11,14,15,16,19,23,25,26,34,36,37,47,48],[2,3,4,5,6,8,9,1,11,15,16,17,24,31],[7,11,15,17,18,27,31,34],[9,13,15,16,12],...
    [4,5,6,8,14,18,19,20,25,26],[2,3,4,5,10,11,17,22,38,47],[4,5,6,7,9,10,13,14,15,18,22],[3,6,7,9,10,12,13,19,21,24,26,29],[6,8,10,13,15,17,28,35,38,45,50],...
    [3,7,11,15,17,22,27,28,40],[6,7,13,14,27,28,29,38,39],[4,5,6,8,10,15,18,23,28,39],[2,3,4,6,7,10,11,13,14,18,20,27,28,30,34],[3,4,6,10,12,14,21,25,26,48],...
    [2,5,6,8,13,15,16,17,21,28,32],[2,5,6,7,8,12,13,22,44],[1,5,7,13,17,18,19,24,32,33,50],[2,7,8,9,10,14,18,18,20,21,22,24],[4,5,6,10,12,13,18,31],...
    [3,4,6,7,8,9,10,11,12,13,14,17,19,20,22,30,35,36,44,61],[4,5,6,7,8,10,12,27,33],[3,6,7,8,12,15,16,19,20,22,25,35],[5,7,10,11,13,14,15,16,18,19,21,22,24,25,28],...
    [5,6,9,11,13],[3,5,6,7,8,9,13,14,16,21,23,28,32,35],[3,4,6,7,8,9,10,11,12,13,16,18,29,39,40,51,59],[4,5,12,13,14,16,20,34,37,40,41],...
    [3,6,9,11,12,16,19,23,24],[3,5,7,8,10,12,14,18,19,20,21,27,28,35],[2,4,5,8,10,13,14,15,16,17]} ;

hfcomps = {[4,8,12,14],[3,7,11,14,25,26],[2,6,8,9,17],[7,11,15],[9,13,18],[6,14,18],[5,11,17],[6,7,10,13],[3,6,7,12,19],[6,10,15,35,45],[3,15,27],[7,27],...
    [4,5,15,16,18],[2,4,6,11,18],[3,6,21],[5,8,16,21],[2,5,7],[5,7,13],[2,7,14],[4,6,10,13],[4,6,9,10,12,13],[4,6,8,10],[3,7,12,15,16],[10,11,14,16],[6,9,11,13],...
    [3,23,28,32],[3,7,8,11,13,16],[4,5,12,14],[6,11,12,19],[7,10,12,19,20],[2,4,8,15]} ; 



for sb=1:length(subs)
    cd(['c:/shared/allres/',subs{sb}]) ; 
    merged = pop_loadset('merged.set') ; 
    comps = load('wcomps.mat') ; 
    comps = comps.wcomps ; weights = comps{1} ; sphere = comps{2} ; 
    acts = weights*sphere*merged.data ; winv = pinv(weights*sphere) ; 
    newmerged = merged ; 
    newmerged.data = acts ; 
    allep = pop_epoch(newmerged,{trigs{1},trigs{3},trigs{4},trigs{5}},[-1,3]) ; 
    for i=1:size(allep.data,1)
        [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',200) ; 
    end
    figure,
    for i=1:64 
        subplottight(10,14,i*2) ; 
        imagesc(squeeze(ersp(i,:,:)),[-4,4]) ; axis xy ; set(gca,'XTick',[],'YTick',[]) ; text(-15,0,num2str(i)) ; 
        subplottight(10,14,i*2-1) ; topoplot(winv(:,i),merged.chanlocs) ; 
    end ; 
    suptitle(subs{sb}) ; 
end
