clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
for sub=1%:length(subs)
cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls
mov = dir('bcgica*set') ; 
currentind = 1 ; 
for m=1:length(mov)
    EEG = pop_loadset(mov(m).name) ; 
    EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
    allsets{m} = EEG ; allinds{m} = currentind:currentind+size(EEG.data,2)-1 ; currentind = currentind + size(EEG.data,2) ; 
    if m==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end
end
mergeica = pop_runica(merged,'runica') ; 


rawts = mean(merged.data([61,29],:),1) - mean(merged.data([62,30],:),1) ;
chan = rawts ; smthchan = smooth(abs(chan),500) ; 

corrs = corr(abs(chan)',abs(mergeica.icaact)') ; 
[sv,si] = sort(corrs,'descend') ; clear allchans 
allchans = mergeica.icaact(si(1:5),:) ; allchans(6,:) = chan ; 

%{
for i=1:64 ; corrs(i) = corr2(smooth(abs(EEG.icaact(i,:)),400),smthchan) ; corrchans(i,:) = smooth(abs(EEG.icaact(i,:)),100) ; end
[sv,si] = sort(corrs,'descend') ; 
sc = mean(corrchans(si(1:3),:)) ; 
plot(sc) ; hline(median(sc)*1.5) ; 
smthchan = sc ; 

plot(smthchan) ; hline(median(smthchan)*2.5) ; 
bads = find(smthchan > median(smthchan)*1.5) ; 
goods = find(smthchan < median(smthchan)*1.5) ; 

bchan = zeros(1,size(EEG.data,2)) ; bchan(bads) = 1 ; bchan([1:round(EEG.srate*20),end-round(EEG.srate*20)]) = 1 ; 
clustbads = bwconncomp(bchan) ; badlist = clustbads.PixelIdxList ; 

gchan = zeros(1,size(EEG.data,2)) ; gchan(bads) = 1 ; 
clustgoods = bwconncomp(gchan==0) ; goodlist = clustgoods.PixelIdxList ; 

for i=1:length(badlist)
    badinds = badlist{i} ;
    badl = length(badinds) ;    
    randind = round(rand*(length(goodlist)-1))+1 ; 
    while length(goodlist{randind}) < badl
        randind = round(rand*(length(goodlist)-1))+1 ; 
    end
    goodinds = goodlist{randind} ; 
    EEG.data(:,badinds) = EEG.data(:,goodinds(1:length(badinds))) ; 
end
%}
end

for cc=1:(size(allchans,1))
    chan = allchans(cc,:) ; disp(cc) ; 
    M = 250 ; Wc = .1 ; k = (1:M-1); s = sin(Wc*k)./ k; c0 = [Wc,s]; A = toeplitz(c0);
    [V,evals] = eig(A); % Only need the principal eigenvector
    [emax,imax] = max(abs(diag(evals)));
    w = ones(1,length(M)) ; 
    g1(:,1:size(w,2)) = w ; 
    halfl = round(M/2) ; wincr = 25 ; icount = 1 ; 
    clear specs  
    for i=halfl+1:wincr:length(chan)-M 
        datai = chan(i-halfl:i+M-halfl-1) ;
        windowed = repmat(datai',[1,size(w,2)]).*w ; 
        f = abs(fft(windowed,[],1)) ; 
        spec = f(1:halfl) ; 
        specs(icount,:) = spec ; 
        icount = icount + 1 ; 
    end
    allspecs(cc,:,:) = specs ; 
end
mspecs = squeeze(mean(allspecs,1)) ; 
plot(mean((mspecs(:,10:end).^2),2)) ; hline(median(mean((mspecs(:,10:end).^2),2))*2)

motion = mean((mspecs(:,10:end).^2),2) ; bw = bwconncomp(imdilate(motion>median(motion)*2.5,strel(ones(5,1)))) ; 
bsizes = cellfun(@length,bw.PixelIdxList) ; 
goods = find(bsizes>10) ; plist = bw.PixelIdxList ; 
zs = zeros(size(motion)) ; for i=1:length(goods) ; zs(plist{goods(i)}) = 1 ; end
bar(zs);  hold on ; plot(mat2gray(motion),'r') ;

reszs = imresize(zs,[size(mergeica.data,2),1]) ; reszs = (reszs~=0) ; 
bar(reszs) ; hold on ; plot(mat2gray(abs(mergeica.icaact(2,:))),'r') ;
figure,subplot(1,2,1) ; plot(rawts((reszs==0))) ; subplot(1,2,2) ; plot(rawts) ; 

nomotion = reszs==0 ; nomotion_clust = bwconncomp(nomotion) ; 
goodlist = nomotion_clust.PixelIdxList ; 

motion = reszs==1 ; motion_clust = bwconncomp(motion) ; 
badlist = motion_clust.PixelIdxList ; 
newEEG = mergeica ; 
for i=1:length(badlist)
    badinds = badlist{i} ;
    badl = length(badinds) ;    
    randind = round(rand*(length(goodlist)-1))+1 ; 
    while length(goodlist{randind}) < badl
        randind = round(rand*(length(goodlist)-1))+1 ; 
    end
    goodinds = goodlist{randind} ; 
    newEEG.data(:,badinds) = newEEG.data(:,goodinds(1:length(badinds))) ; 
end
newmergeica = pop_runica(newEEG,'runica') ; 
plot(mergeica.data(45,:),'r') ; hold on ; plot(newEEG.data(45,:)) ; 

for i=1:length(allsets)
    eegi = allsets{i} ; size(eegi.data), size(allinds{i})
    eegi.data = newmergeica.data(:,allinds{i}) ;  
    eegi = ica_applyweights(eegi,newmergeica) ; 
    pop_saveset(eegi,strrep(mov(i).name,'bcgica','mergeica')) ; 
end

