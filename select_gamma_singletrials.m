clear all ; close all ; 
cd C:\shared\allres
subs = {'alex','alexandra3','audrey','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','leila','lisa','marc',...
    'marie','mathieu','maxime','menglu','mingham','olga','patricia','po','russell','suhan2','sunachakan','tah','tegan2','vincent'} ; 
trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'} ; 

hfcomps = {[4,8,12,14],[3,7,11,14,25,26],[2,6,8,9,17],[7,11,15],[9,13,18],[6,14,18],[5,11,17],[6,7,10,13],[3,6,7,12,19],[6,10,15,35,45],[3,15,27],[7,27],...
    [4,5,15,16,18],[2,4,6,11,18],[3,6,21],[5,8,16,21],[2,5,7],[5,7,13],[2,7,14],[4,6,10,13],[4,6,9,10,12,13],[4,6,8,10],[3,7,12,15,16],[10,11,14,16],[6,9,11,13],...
    [3,23,28,32],[3,7,8,11,13,16],[4,5,12,14],[6,11,12,19],[7,10,12,19,20],[2,4,8,15]} ; 

for sb=1:length(subs)
    cd(['c:/shared/allres/',subs{sb}]) ; 
    merged = pop_loadset('merged.set') ; 
    merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61) ; 
    comps = load('wcomps.mat') ; 
    comps = comps.wcomps ; weights = comps{1} ; sphere = comps{2} ; winv = pinv(weights*sphere) ; 
    acts = weights*sphere*merged.data ; 
    zs = zeros(1,64) ; zs(hfcomps{sb}) = 1 ; bads = find(zs==0) ; 
    %acts(bads,:) = 0 ; 
    merged.icaact = acts ; merged.icaweights = weights ; merged.icasphere = sphere ; merged.icawinv = winv ; merged.icachansind = 1:64 ; 
    merged = pop_subcomp(merged,bads) ; 
    
    freqs = 1:4:120 ; 
    for tr=1:length(trigs)
        allep = pop_epoch(merged,{trigs{tr}},[-1.5,3]) ; 
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
        allfreqeps(sb,tr,:,:,:) = mresfreqeps ; 
    end
    %for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mresfreqeps(i,:,:)),[-.2,.2]) ; title(i) ; end
    %figure,topoplot(squeeze(mean(mean(mresfreqeps(:,freqs>50 & freqs<80,allep.times>0 & allep.times<2000),2),3)),merged.chanlocs,'electrodes','numbers')  ;title(subs{sb}) ; 
    %alltopos(sb,:) = squeeze(mean(mean(mresfreqeps(:,freqs>50 & freqs<80,allep.times>0 & allep.times<2000),2),3)) ;
    figure,topoplot(squeeze(mean(mean(mean(allfreqeps(sb,:,:,freqs>50 & freqs<90,allep.times>0 & allep.times<2000),4),5),2)),merged.chanlocs) ; 
    
end

for sb=1:length(subs)
    cd(['c:/shared/allres/',subs{sb}]) ; 
    afreps = squeeze(allfreqeps(sb,:,:,:,:)) ; 
    save('afreps','afreps') ; 
end
