clear all ; close all 
meanmotion = zeros(7,6) ; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie'} ; 
for s=1:length(subs)
cd(['c:/shared/badger_eeg/',subs{s}]) ; 
grads = dir('grad_*set') ; 
for g=1:length(grads)
EEG = pop_loadset(grads(g).name) ; 
filtdat = eegfiltfft(EEG.data,EEG.srate,80,128) ; 
ts = (smooth(std((filtdat([1:31,33:end],:)),0,1),100)) ;    
tslp = eegfiltfft(ts',256,0.01,100) ; 

% perform iterative z-score thresholding:
tempts = tslp ; zts = zeros(1,length(tslp)) ; 
for i=1:50
   zs = zscore(tempts) ; 
   [sv,si] = sort(zs,'descend') ; 
   tempts(si(1:100)) = 0 ; 
   zts(si((sv>3))) = 1 ; 
end
figure, 

% find and eliminate small z clusters (<40 continuous time points)
bwts = bwconncomp(zts>0) ; 
inds = bwts.PixelIdxList ; 
indlengths = cellfun(@length,inds) ; 
goodinds = find(indlengths>25) ; 
findinds = 1:length(zts) ; % to find the indices (no need for bounds-checking) 
srate = EEG.srate ; 
newzs = zeros(1,length(zts)) ; 
for i=1:length(goodinds)
    newzs(bwts.PixelIdxList{goodinds(i)}) = 1 ; 
    minind = round(min(bwts.PixelIdxList{goodinds(i)})) ;
    maxind = round(max(bwts.PixelIdxList{goodinds(i)})) ; 
    surroundinds = find(findinds < maxind+srate/2 & findinds > minind-srate/2) ; 
    newzs(surroundinds) = 1 ; 
end

% get all the motion epochs
clustinds = bwconncomp(newzs) ; 
for i=1:length(clustinds.PixelIdxList)
    epochs{i} = tslp(clustinds.PixelIdxList{i}) ; 
end

% get a random sample of data epochs
nsamples = 500 ; 
maxsamplelength = max(cellfun(@length,epochs)) ; 
clear samples 
for i=1:nsamples
    rndind = round(rand(1)*(length(zts)-maxsamplelength-1)) ; 
    samples(i,:) = tslp(rndind:rndind+maxsamplelength) ; 
end

% fill up the number line with histogram indices
sumsamples = sum(diff(samples,1,2).^2,2) ; 
badsamps = find(zscore(log(sumsamples))>0) ; 
samples(badsamps,:) = [] ;
diffsamples = diff(samples,1,2) ; 

% find a way to automatically get these parameters from the data (zscore?)



for i=1:length(epochs)
    epochi = epochs{i} ; 
    diffepochi = diff(epochi,1,2) ; 
    samphist = histc(diffsamples(:,1:length(diffepochi)),[-.03:.001:.03],2) ;
    cmat = corr(samphist') ;
    meancorrs = mean(cmat,2) ;   
    histi = histc(diffepochi,[-.03:.001:.03],2) ;
    corrsi = corr(histi',samphist') ; 
    %plot(meancorrs) ; hold on ; plot(corrsi,'r') ; 
    %title(num2str(mean(meancorrs-corrsi')./(std(meancorrs) + std(corrsi)))) ; 
    allshists(i,:) = histi ; 
    ttests(i) = mean(meancorrs-corrsi')./(std(meancorrs) + std(corrsi)) ; 
    
end

% remove the false positives
clustlist = clustinds.PixelIdxList ; 
nofalsepos = zeros(1,length(newzs)) ; 
for i=1:length(clustlist) 
    if ttests(i) > 2
        nofalsepos(clustlist{i}) = 1 ; 
    end
end

%for i=1:length(epochs) ; plot(epochs{i}) ; hold on ; end
%for i=1:length(epochs) ; subplot(ceil(sqrt(length(epochs))),ceil(sqrt(length(epochs))),i) ; plot(epochs{i}) ; title(ttests(i)) ; end
%bar(newzs) ; 
%bar(newzs,'r') ;hold on ;  plot(mat2gray(tslp),'k'); ylim([0,.15])
meanmotion(s,g) = sum(newzs) ; 
end
end
