clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
for sub=3%:length(subs)
cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls
mov = dir('weights*movie*set') ; 
EEG = pop_loadset(mov(1).name) ; 
EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;

rawts = mean(EEG.data([61,29],:),1) - mean(EEG.data([62,30],:),1) ;
chan = rawts ; smthchan = smooth(abs(chan),500) ; 

for i=1:64 ; corrs(i) = corr2(smooth(abs(EEG.icaact(i,:)),400),smthchan) ; corrchans(i,:) = smooth(abs(EEG.icaact(i,:)),100) ; end
[sv,si] = sort(corrs,'descend') ; 
sc = mean(corrchans(si(1:3),:)) ; 
plot(sc) ; hline(median(sc)*1.5) ; 
smthchan = sc ; 

plot(smthchan) ; hline(median(smthchan)*2.5) ; 
bads = find(smthchan > median(smthchan)*2) ; 
goods = find(smthchan < median(smthchan)*2) ; 

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
EEG2 = EEG ; 
EEG2 = pop_runica(EEG2,'runica') ; 
tp(EEG2) ; 
goods = [4,5,6,7,9,11,12,14,16,19,21,22,23,32,30,43,45] ; zs = zeros(1,64) ; zs(goods) = 1 ; bads = find(zs==0) ; 
subeeg = pop_subcomp(EEG2,bads)  ;

end


M = 250 ; Wc = .1 ; k = (1:M-1); s = sin(Wc*k)./ k; c0 = [Wc,s]; A = toeplitz(c0);
[V,evals] = eig(A); % Only need the principal eigenvector
[emax,imax] = max(abs(diag(evals)));
w = ones(1,length(M)) ; 
g1(:,1:size(w,2)) = w ; 
halfl = round(M/2) ; wincr = 25 ; 
chans = [1:64] ; clear specs allchans 
for c=1:length(chans) ; disp(['chan=',num2str(c)]) ; 
chan = subeeg.data(chans(c),:) ; icount = 1 ; 
for i=halfl+1:wincr:length(chan)-M 
    datai = chan(i-halfl:i+M-halfl-1) ;
    windowed = repmat(datai',[1,size(w,2)]).*w ; 
    f = abs(fft(windowed,[],1)) ; 
    spec = f(1:halfl) ; 
    specs(icount,:) = spec ; 
    icount = icount + 1 ; 
end
allchans(c,:,:) = specs ; 
end

cd(['c:/shared/newbadger_mri/',subs{sub}])
corrs = load_untouch_nii('corrs.nii.gz') ; 
rest = load_untouch_nii('bp_reg_topup_mc_retino_movie.nii.gz') ; 
corrvox = find(medfilt3(corrs.img)>.125) ; [i1,i2,i3] = ind2sub(size(corrs.img),corrvox) ; ts = zeros(length(i1),size(rest.img,4)) ; 
for i=1:length(i1) ; ts(i,:) = squeeze(rest.img(i1(i),i2(i),i3(i),:)) ; end
mts = mean(ts(:,1:end-1)) ; 

lats = cell2mat({EEG.urevent.latency}) ; 
labs = {EEG.urevent.type} ; 
r128s = find(strcmp('R128',labs)) ; 
r128lats = round(lats(r128s)/wincr) ; 

tracts = allchans(:,r128lats(1):r128lats(end),:) ; 
reshape_sz = length(r128s)-1 ; restracts = zeros(size(tracts,1),reshape_sz,size(tracts,3)) ; 
for i=1:size(tracts,1) ; restracts(i,:,:) = imresize(medfilt2(squeeze(((tracts(i,:,:)))),[5,5]),[reshape_sz,size(tracts,3)]) ; end
filteeg = zeros(size(restracts)) ; 
for i=1:64 ; rts = eegfiltfft(squeeze(restracts(i,:,:))',1/0.693,0.02,1.5) ; filteeg(i,:,:) = rts' ; end

wincr = 25 ; wsize = 150 ; tcount = 1 ; clear xcorrs ; 
for t=35:wincr:size(mts,2)-(wsize+35)
    for i=1:64 ; 
        for j=1:125
            fmri = mts(t:t+wsize) ;
            eeg = squeeze(restracts(i,t:t+wsize,j)) ; 
            xcorrs(tcount,i,j,:) = xcorr(fmri,eeg,25,'coeff') ; 
        end
    end 
    tcount = tcount + 1 ; 
end

tp(EEG) ; 
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(xcorrs(:,i,:,:),1)),[-.2,.2]) ; title(i) ; axis xy ; end
chans = [23,64,37,19] ;
figure,imagesc(squeeze(mean(mean(xcorrs(:,chans,:,:),1),2))) ; axis xy ; 

figure,topoplot(squeeze(mean(xcorrs(:,:,10,36),1)),EEG.chanlocs,'electrodes','numbers') ; 






