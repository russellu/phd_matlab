clear all ; close all
cd C:\shared\allres
subs = {'alex','alexandra3','audrey','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','leila','lisa','marc',...
    'marie','mathieu','maxime','menglu','mingham','patricia','po','russell','suhan2','sunachakan','tah','tegan2','vincent'} ; 
badchans2 = {[22,17],[28,12,16],[42,22],[53,26],[17,45,63],[16,12],[17,22],[28,22],[32,22,46],[22,33,41],[3,22,44],[12,46],[22],[32,28,12,2],[22,32,28],[17,35,34],[29,32,60],[22,55,42,12],[32,28],[22,53,48,42],[17],[17,2],[34,63,49],[22,28,2],[17,28],[32,17],[22,32,19,14,28],[32],[22,28,16,17],[32,25,64]} ; 
swaps = [1,1,0,0,0,0,0,1,1,1,1,0,0,1,1,1,0,0,1,1,0,0,0,1,0,0,1,0,0,1] ; 
trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'} ; 

for sb=1:length(subs)
    cd(['c:/shared/allres/',subs{sb}]) ; 
    ls 
    merged = pop_loadset('cleanfilt.set') ; 
    %merged = pop_resample(merged,256) ; 
    %pop_saveset(merged,'cleanfilt.set') ; 
    
    fullica = load('fullica.mat') ; fullica = fullica.fullica ; weights = fullica{1} ; sphere = fullica{2} ; 
    for tr=1:length(trigs)
        epica = pop_epoch(merged,{trigs{tr}},[-1,3]) ; clear ersp ; 
        for i=1:64 ; disp(i)
            for j=1:size(epica.data,3)
                [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(i,:,j)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',NaN,'verbose','off','timesout',25) ; 
            end
        end
        allersp{tr} = ersp ; 
    end
    
    for i=1:length(allersp)
        ersp = allersp{i} ; 
        zs = (zscore(sum(squeeze(mean(std(ersp,0,4),3))),1)) ; 
        [~,si{i}] = sort(zs,'ascend') ; 
    end
    save('si','si') ; 
    
    merged.data = weights*sphere*merged.data ; 
    clear bersp ; 
    for tr=1:length(trigs)
        ep = pop_epoch(merged,{trigs{tr}},[-1,3]) ; 
    for i=1:64 ; 
        [bersp(tr,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.data(i,:,si{tr}(1:end-10))),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                            'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',100) ; 
    end
    end
    save('bersp','bersp') ; 
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(bersp(:,i,:,:),1)),[-2,2]) ; end
    
    clear ersp ; 
    for tr=1:length(trigs)
        ep = pop_epoch(merged,{trigs{tr}},[-1,3]) ; 
    for i=1:64 ; 
        [ersp(tr,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.data(i,:,si{tr}(1:end-10))),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                            'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off','timesout',100) ; 
    end
    end
    save('ersp','ersp') ; 
    save('times','times') ; save('freqs','freqs') ; 

end