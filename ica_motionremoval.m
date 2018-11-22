clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
for sub=9%:length(subs)
cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls
mov = dir('bcgica*rest*set') ; 
EEG = pop_loadset(mov(1).name) ; 
EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
rawts = mean(EEG.data([61,29],:),1) - mean(EEG.data([62,30],:),1) ;
smth = smooth(rawts.^2,300) ; 
ica = pop_runica(EEG,'runica','maxsteps',128) ; 
end

M = 250 ; Wc = .1 ; k = (1:M-1); s = sin(Wc*k)./ k; c0 = [Wc,s]; A = toeplitz(c0);
[V,evals] = eig(A); % Only need the principal eigenvector
[emax,imax] = max(abs(diag(evals)));
%w = V(:,end-2:end) ;clear g1 
w = ones(1,length(M)) ; 
g1(:,1:size(w,2)) = w ; 
ica.icaact(64,:) = rawts ; 
halfl = round(M/2) ; wincr = 25 ; 
chans = [1:64] ; clear specs allchans 
for c=1:length(chans) ; disp(['chan=',num2str(c)]) ; 
chan = ica.icaact(chans(c),:) ; icount = 1 ; 
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
basechans = allchans - repmat(mean(allchans,2),[1,size(allchans,2),1]) ; 
%{
resbase = reshape(basechans,[size(basechans,1),numel(basechans(1,:,:))]) ; 
[weights,sphere] = runica(resbase(:,1:3:end),'maxsteps',128) ; 
acts = weights*sphere*resbase ; resacts = reshape(acts,size(basechans)) ; 
%}
clear allds ; 
for i=1:125
[s,v,d] = svd(allchans(:,:,i)) ;
s(:,1:5) = 0 ; 
allds(i,:,:) = d(:,1:5) ; 
reschan = s*v*d' ; 
reschans(:,:,i) = reschan ; disp(i) ; 
end




%{
pt = (squeeze(mean(mean((allchans),1),3))') + squeeze(mean(allchans(64,:,:),3))'  ; 
mchans = squeeze(mean(allchans,3)) ; 
corrs = corr(mchans') ; 
[sv,si] = sort(corrs(64,:).^2,'descend') ; 
meanmotion = mean(squeeze(mean(allchans(si(1:5),:,:),1))') ; 
meanfreqs = squeeze(mean(allchans(si(1:5),:,:),1)) ; 
badmotion = (meanmotion > median(meanmotion)*1.5) ; 
badcomp = bwconncomp(badmotion) ; badlist = badcomp.PixelIdxList ; 
badszs = cellfun(@length,badlist) ; goodclusts = badlist(badszs>5) ; 
zs = zeros(size(meanmotion)) ; for i=1:length(goodclusts) ; zs(goodclusts{i}) = 1 ; end
negclust = bwconncomp(zs==0) ; negszs = cellfun(@length,negclust.PixelIdxList) ; 
toosmall = find(negszs<=7) ; for i=1:length(toosmall) ; zs(negclust.PixelIdxList{toosmall(i)}) = 1 ; end
zs = imdilate(zs,strel(ones(1,3))) ; 
figure,bar(zs) ; hold on ; plot(mat2gray(meanmotion),'r') ;
interpchans = allchans ; 
for i=1:64 ; disp(i) ; 
meanfreqs = double(squeeze(allchans(i,:,:))) ; 
[gx,gy] = meshgrid(1:size(meanfreqs,2),1:size(meanfreqs,1)) ;
goodinds = find(zs==0) ; badinds = find(zs==1) ;
goodx = gx(goodinds,:) ; goody = gy(goodinds,:) ; goodv = meanfreqs(goodinds,:) ; 
badx = gx(badinds,:) ; bady = gy(badinds,:) ; 
v1 = griddata(goodx,goody,goodv,badx,bady) ; 
meanfreqs(badinds,:) = v1 ; 
interpchans(i,:,:) = meanfreqs ; 
end
%}

cd(['c:/shared/newbadger_mri/',subs{sub}])
corrs = load_untouch_nii('corrs.nii.gz') ; 
rest = load_untouch_nii('bp_reg_topup_mc_retino_rest.nii.gz') ; 
corrvox = find(medfilt3(corrs.img)>.125) ; [i1,i2,i3] = ind2sub(size(corrs.img),corrvox) ; ts = zeros(length(i1),size(rest.img,4)) ; 
for i=1:length(i1) ; ts(i,:) = squeeze(rest.img(i1(i),i2(i),i3(i),:)) ; end
mts = mean(ts(:,1:end-1)) ; 

lats = cell2mat({EEG.urevent.latency}) ; 
labs = {EEG.urevent.type} ; 
r128s = find(strcmp('R128',labs)) ; 
r128lats = round(lats(r128s)/wincr) ; 

tracts = reschans(:,r128lats(1):r128lats(end),:) ; 
reshape_sz = length(r128s)-1 ; restracts = zeros(size(tracts,1),reshape_sz,size(tracts,3)) ; 
for i=1:size(tracts,1) ; restracts(i,:,:) = imresize(medfilt2(squeeze(((tracts(i,:,:)))),[5,5]),[reshape_sz,size(tracts,3)]) ; end
filteeg = zeros(size(restracts)) ; 
for i=1:64 ; rts = eegfiltfft(squeeze(restracts(i,:,:))',1/0.693,0.02,1.5) ; filteeg(i,:,:) = rts' ; end


meanres = squeeze(mean(restracts(chans,:,:))) ; 
wincr = 10 ; wsize = 50 ; tcount = 1 ; clear xcorrs ; 
for t=35:wincr:size(mts,2)-(wsize+35)
    for i=1:64 ; 
        for j=1:125
            fmri = mts(t:t+wsize) ;
            eeg = squeeze(restracts(i,t:t+wsize,j)) ; 
            %eeg = meanres(t:t+wsize,j) ; 
            xcorrs(tcount,i,j,:) = xcorr(fmri,eeg,20,'coeff') ; 
        end
    end 
    tcount = tcount + 1 ; 
end
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(xcorrs(:,i,:,:),1)),[-.25,.25]) ; title(i) ; axis xy ; end

%{
chans = [5,7,8,14] ;
figure,imagesc(squeeze(mean(mean(xcorrs(:,chans,:,:),1),2)),[-.25,.25]) ;  axis xy ; 


figure,for i=1:size(xcorrs,1) ; subplot(6,11,i) ; imagesc(squeeze(mean(xcorrs(i,chans,:,:)))) ; axis xy ; title(i)  ;end ; 

subplot(2,2,1) ;imagesc(squeeze(mean(mean(xcorrs(:,chans,:,:),1),2)),[-.3,.3]) ;  axis xy ; vline(20,'k') ; 
subplot(2,2,2) ;
plot(squeeze(mean(mean(mean(xcorrs(:,chans,12:18,:),1),2),3))) ;vline(20,'k') ; hline(0,'k') ; xlim([0,41])

%}

