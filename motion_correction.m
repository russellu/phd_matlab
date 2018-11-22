clear all ; close all 
meanmotion = zeros(7,6) ; 
subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 
for s=1:length(subs) ; 
cd(['c:/shared/badger_eeg/',subs{s}]) ; 
grads = dir('*rest*Pulse*set') ; 
EEG = pop_loadset(grads(1).name) ; 
edat = EEG.data([1:31,33:64],:) ; edat(:,[1:10000,size(edat,2)-10000:size(edat,2)]) = 0 ; 
filtdat = eegfiltfft(EEG.data,EEG.srate,30,120) ; 
[weights,sphere] = runica(filtdat(:,10000:end-10000),'maxsteps',150) ;
acts = weights*sphere*filtdat ; 
icastd = (smooth(std(acts(:,1:end),0,1),100)) ; 
zthresh = icastd ; bads = zeros(1,length(zthresh)) ; 
for i=1:50
    zinds = find(zscore(zthresh)>2) ; 
    bads(zinds) = 1 ; 
    zthresh(zinds) = 0 ; 
end
dilbads = (imdilate(bads,strel(ones(1,round(EEG.srate)*2)))) ; 
figure,subplot(1,2,1) ; bar(dilbads) ; hold on ; plot(mat2gray(icastd),'r') ;
clustbads = bwconncomp(dilbads) ; 
for i=1:length(clustbads.PixelIdxList) ; 
 %  plot(icastd(clustbads.PixelIdxList{i})) ; hold on ;  
   badepochs{i} = icastd(clustbads.PixelIdxList{i}) ; 
   bars(i) = sum(zscore(badepochs{i})>1)./length(badepochs{i}) ; 
end
EEG.clustlist = clustbads.PixelIdxList ; 
badts = zeros(1,length(dilbads)) ; 
for i=1:length(clustbads.PixelIdxList) ; badts(clustbads.PixelIdxList{i}) = 1 ; end ; 
badts = find(badts == 1) ; 
EEG.badts = badts ; 
pop_saveset(EEG,grads(1).name) ; 
%{
for i=1:length(badepochs)
    subplot(ceil(sqrt(length(badepochs))),ceil(sqrt(length(badepochs))),i) ; 
    plot(badepochs{i}) ; title(bars(i)) ; 
end

newbads = zeros(1,length(dilbads)) ; 
for i=1:length(clustbads.PixelIdxList)
    if bars(i) > 2.5
        newbads(clustbads.PixelIdxList{i}) = 1 ; 
    end
end
subplot(1,2,2) ; bar(newbads) ; hold on ; plot(mat2gray(icastd),'r') ;
%}
end
    