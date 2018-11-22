clear all ; close all
cd C:\shared\allres
subs = {'alex','alexandra3','audrey','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','leila','lisa','marc',...
    'marie','mathieu','maxime','menglu','mingham','olga','patricia','po','russell','suhan2','sunachakan','tah','tegan2','vincent'} ; 
badchans = {[49,54],[35,44,48],[22,42],[42],[7,17,45],[12,16,21],[17],[54,60],[54,64],[1,54],[12,35,54],[42],[22],[60,64],[54,64],[8,49],[45],[22],[13,43],...
    [16,21,54],[17],[42,45],[10,11,35,50],[34,63],[54],[17],[17,32],[51,54,60,64],[32],[8,22,42],[64]} ;

badchans2 = {[22,17],[28,12,16],[42,22],[53,26],[17,45,63],[16,12],[17,22],[28,22],[32,22,46],[22,33,41],[3,22,44],[12,46],[22],[32,28,12,2],[22,32,28],[17,35,34],[29,32,60],[22,55,42,12],[32,28],[22,53,48,42],[17],[17,2],[17,16,22],[34,63,49],[22,28,2],[17,28],[32,17],[22,32,19,14,28],[32],[22,28,16,17],[32,25,64]} ; 
swaps = [1,1,0,0,0,0,0,1,1,1,1,0,0,1,1,1,0,0,1,1,0,0,1,0,1,0,0,1,0,0,1] ; 

for sb=1:length(subs)
    cd(['c:/shared/allres/',subs{sb}]) ; 
    ls 
    resamps = dir('resamp*set') ; 
    for res=1:length(resamps)
        eeg = pop_loadset(resamps(res).name) ;
        if swaps(sb) == 1 ; newdat = eeg.data(1:32,:) ; eeg.data(1:32,:) = eeg.data(33:64,:) ; eeg.data(33:64,:) = newdat ; end 
        zs(res,:) = zscore(sum(abs(diff(eeg.data,1,2)),2)) ; 
        if res==1 ; merged = eeg ; else merged = pop_mergeset(eeg,merged) ; end
    end
    [sv,si] = sort(max(zs,[],1),'descend') ; 
    allzs(sb,:) = max(zs,[],1) ; 
    %figure('units','normalized','outerposition',[0 0 1 1]),subplot(1,2,1) ; bar(max(zs,[],1)) ; suptitle([subs{sb},',',num2str(si(1)),',',num2str(si(2)),',',num2str(si(3)),',',...
    %    num2str(si(4)),',',num2str(si(5)),',',num2str(si(6)),',',num2str(si(7)),',',num2str(si(8)),',',num2str(si(9)),',',num2str(si(10))]) ; 
    %subplot(1,2,2) ; topoplot(max(zs,[],1),merged.chanlocs,'maplimits',[-5,5],'electrodes','numbers') ; 
end



