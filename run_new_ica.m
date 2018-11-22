clear all ; close all ; 
cd C:\shared\allres
subs = {'alex','alexandra3','audrey','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','leila','lisa','marc',...
    'marie','mathieu','maxime','menglu','mingham','olga','patricia','po','russell','suhan2','sunachakan','tah','tegan2','vincent'} ; 
trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'} ; 


for sb=22%:length(subs)
    cd(['c:/shared/allres/',subs{sb}]) ; 
    merged = pop_loadset('merged.set') ; 
    merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61) ; 
    newmerged = merged ; newmerged.data = eegfiltfft(merged.data,merged.srate,40,90) ; 
    newmerged = pop_epoch(newmerged,trigs,[-1,3]) ; 
    newmerged = pop_runica(newmerged,'runica','maxsteps',128) ; 
    rawmerged = merged ; 
    merged.data = newmerged.icaweights*newmerged.icasphere*merged.data ; 
    allep = pop_epoch(merged,{trigs{1},trigs{5}},[-1,3]) ; 
    for i=1:size(allep.data,1)
        [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',200) ; 
    end
    figure,
    for i=1:64 
        subplottight(10,14,i*2) ; 
        imagesc(squeeze(ersp(i,:,:)),[-4,4]) ; axis xy ; set(gca,'XTick',[],'YTick',[]) ; text(-15,0,num2str(i)) ; 
        subplottight(10,14,i*2-1) ; topoplot(newmerged.icawinv(:,i),merged.chanlocs) ; 
    end ; 
    
    
    
    goodcs = [20,30,35] ; chans = zeros(1,64) ; chans(goodcs) = 1 ; bads = find(chans==0) ; 
    rawmerged.icaact = newmerged.icaweights*newmerged.icasphere*rawmerged.data ; 
    rawmerged.icachansind = newmerged.icachansind ; rawmerged.icaweights = newmerged.icaweights ; rawmerged.icasphere = newmerged.icasphere ; 
    rawmerged.icawinv = pinv(rawmerged.icaweights*rawmerged.icasphere) ; 
    rawmerged = pop_subcomp(rawmerged,bads) ; 
    %{
    testelecs = 1:10:64 ; 
    clear allersp ; 
    for tr=1:length(trigs)
        allep = pop_epoch(subbed,{trigs{tr}},[-1,3]) ; 
        ersp = zeros(length(testelecs),size(allep.data,3),60,100) ; 
        for i=1:length(testelecs) ; disp(i) ; 
            for j=1:135
            [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(testelecs(i),:,j)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off','timesout',100) ; 
            end
        end
        allersp{tr} = ersp ; 
    end
    for i=1:length(allersp)
        ersp = allersp{i} ; 
        %subplot(2,3,i) ; bar(zscore(squeeze(mean(std(mean(ersp,3),0,4),1))))
        bads = find(zscore(squeeze(mean(std(mean(ersp,3),0,4),1)))>1.5) ; 
        goods = zeros(1,size(allersp{i},2)) ; goods(bads) = 1 ; goodtrials{i} = find(goods==0) ; 
    end
    %}
    %{
    clear ersp ; 
    for tr=1:length(trigs)
        allep = pop_epoch(subbed,{trigs{tr}},[-1,3]) ; 
        for i=1:64 ; 
            [ersp(tr,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,goodtrials{tr})),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',100) ; 
        end
    end
    elecs = [29,30,31,60,61,62,63,64,23,56,24,57,25,58,26,59,27] ; 
    figure,for i=1:6 ; subplot(2,3,i) ; imagesc(times,freqs,squeeze(mean(ersp(i,elecs(1:3),:,:),2)),[-5,5]) ; vline([0,2],'k') ; axis xy ; end  ;suptitle(subs{sb}) ; 
    %}    
    freqs = 1:2:120 ;  
    for tr=1%:length(trigs)
        allep = pop_epoch(rawmerged,{trigs{tr}},[-1.5,3]) ; 
        resep = reshape(allep.data,[64,size(allep.data,2)*size(allep.data,3)]) ; 
        freqeps = zeros(size(resep,1),length(freqs),size(resep,2)) ; 
        for f=1:length(freqs)
            freqeps(:,f,:) = eegfiltfft(resep,merged.srate,freqs(f)-2,freqs(f)+2) ; 
            
        end
        resfreqeps = zeros(size(allep.data,1),size(freqeps,2),size(allep.data,2),size(allep.data,3)) ; 
        for i=1:size(freqeps,2)
            resfreqeps(:,i,:,:) = reshape(squeeze(freqeps(:,i,:)),size(allep.data)) ; 
        end
        clear freqeps ; 
        resfreqeps = abs(resfreqeps) ; 
        resfreqeps = resfreqeps - repmat(mean(resfreqeps(:,:,allep.times<0 & allep.times>-800,:),3),[1,1,size(resfreqeps,3),1]) ; 
        stdfreqeps = squeeze(std(mean(resfreqeps(:,20:end,:,:),2),0,3)) ; 
        [sv,si] = sort(mean(stdfreqeps),'descend') ; 
        mresfreqeps = squeeze(mean(resfreqeps(:,:,:,si(15:end)),4)) ; 
    end
    for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mresfreqeps(i,:,:)),[-.2,.2]) ; title(i) ; end
    figure,topoplot(squeeze(mean(mean(mresfreqeps(:,30:40,allep.times>0 & allep.times<2000),2),3)),merged.chanlocs,'electrodes','numbers','maplimits',[-.02,.02])  ;title(subs{sb}) ; 

end

