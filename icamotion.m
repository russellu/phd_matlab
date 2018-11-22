clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
for sub=7%:length(subs)
cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls
mov = dir('bcgica*movie*set') ; 
EEG = pop_loadset(mov(1).name) ; 
EEG.data = eegfiltfft(EEG.data,250,1,250) ;
EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;

%{
newz = EEG.data ; eeg4 = EEG ; 
rawts = mean(newz([61,29],:),1) - mean(newz([62,30],:),1) ;
chan = rawts ;
wsize = 200 ; lchan = length(chan) ;
spec = zeros(length(chan),100) ; 
for i=1:size(newz,2)-wsize ; disp(i/lchan) ; 
    fi = abs(fft(chan(i:i+wsize))) ; 
    spec(i,:) = fi(1:100) ;
end
smth = smooth(mean(spec(:,70:100),2).^2,350) ; 
figure,plot(smth) ; hline(median(smth)*3)
%}

end
%{
for i=1%:length(allspecs) ; 
    %subplot(3,3,i) ;
    figure,imagesc(spec') ; axis xy ;    
    filti = smooth(mean((spec(:,1:45)),2),150) ;
    %plot(filti) ; hline(median(filti)*2.5,'r') ;
    numbads(i) = numel(find(filti>median(filti)*2.5)) ; 
    bin = (filti> median(filti)*2);  
    dilfilt = imdilate(bin,strel(ones([950,1]))) ; 
    bw = bwconncomp(dilfilt==0) ; 
    plist = bw.PixelIdxList ; 
    for j=1:length(plist)
       listj = plist{j} ; 
       if length(listj) < 250
          dilfilt(listj) = 1 ;  
       end   
    end
    figure,bar(dilfilt*5000) ; hold on ; plot(filti,'r') ; 
end
%}


eeg2 = EEG ; eeg2.data = eegfiltfft(EEG.data,250,1,125) ;%+ eegfiltfft(eeg4.data,250,10,14) ; 
ica2 = pop_runica(eeg2,'runica') ; 
ica2 = ica_applyweights(EEG,ica2) ; tp(ica2) ;

M = 250 ; Wc = .1 ; k = (1:M-1); s = sin(Wc*k)./ k; c0 = [Wc,s]; A = toeplitz(c0);
[V,evals] = eig(A); % Only need the principal eigenvector
[emax,imax] = max(abs(diag(evals)));
w = V(:,end-5:end) ;clear g1 
w = ones(1,length(M)) ; 
g1(:,1:size(w,2)) = w ; 
halfl = round(M/2) ; 
chans = [1:64] ; clear specs allchans 
for c=1:length(chans) ; disp(['chan=',num2str(c)]) ; 
chan = ica2.icaact(chans(c),:) ; icount = 1 ; 
for i=1:10:length(chan)-M 
    datai = chan(i:i+M-1) ;
    windowed = repmat(datai',[1,size(w,2)]).*w ; 
    f = abs(fft(windowed,[],1)) ; 
    spec = f(1:halfl) ; 
    specs(icount,:) = spec ; 
    icount = icount + 1 ; 
end
allchans(c,:,:) = specs ; 
end

goodcomps = [9,20,21] ; md = log(squeeze(mean(allchans(goodcomps,:,:),1))) ; 
md = md - repmat(mean(md,1),[size(md,1),1]) ; 
md = medfilt2(md,[5,5]) ; 

evts = {EEG.urevent.type} ; 
lats = cell2mat({EEG.urevent.latency}) ; 
s1s = find(strcmp('S  1',evts)) ; 
r128s = find(strcmp('R128',evts)) ; 
s1lats = lats(s1s) ; r128lats = lats(r128s) ; 
windowlats = s1lats/10 ; windowr128s = r128lats/10 ; 
ntrs = length(r128s)-1 ; 

st = windowr128s(1) ; en = windowr128s(end) ; 
newmd = md(st:en,:)' ; motion_resmd = imresize(newmd,[halfl,ntrs]) ; 



resbase = reshape(allchans,[64,size(allchans,2)*size(allchans,3)]) ; 
[weights,sphere] = runica(resbase(:,1:5:end),'maxsteps',128) ; 

% remove the motion components
acts = weights*sphere*resbase ; 
resacts = reshape(acts,size(allchans)) ; 
bads = zeros(1,64) ; 
for i=1:64 ; 
    vec = (squeeze(mean(resacts(i,:,:),3))') ; 
    if skewness(vec) > 4
        bads(i) = 1 ; 
    end
end
bads = find(bads==1) ;
acts(bads,:) = 0 ; 
winv = pinv(weights*sphere) ; 
invacts = winv*acts ; 
resinvacts = reshape(invacts,size(allchans)) ; resinvacts = resinvacts - min(resinvacts(:)) ; 
md = log(squeeze(mean(resinvacts(goodcomps,:,:),1))) ; 
md = md - repmat(mean(md,1),[size(md,1),1]) ; 
md = medfilt2(md,[5,5]) ; 
newmd = md(st:en,:)' ; resmd = imresize(newmd,[halfl,ntrs]) ; 
resmd = eegfiltfft(resmd,1/0.693,0.02,1.5) ; 

cd(['c:/shared/newbadger_mri/',subs{sub}]) ; ls
bp = load_untouch_nii('bp_reg_topup_mc_retino_movie.nii.gz') ; bpimg = bp.img ; 

corr = load_untouch_nii('corrs.nii.gz') ; 
mask = medfilt3(corr.img > .25) ; inds = find(mask==1) ; [i1,i2,i3] = ind2sub(size(mask),inds) ; 
for i=1:length(inds)
    ts(i,:) = squeeze(bpimg(i1(i),i2(i),i3(i),:)) ; 
end
mts = mean(ts,1) ; mts = mts(1:ntrs) ; 

%{
clear xcorrs ; 
corrwsize = 120 ; 
corrincr = 10 ; count = 1 ; 
for c=35:corrincr:size(resmd,2)-(corrwsize+35)
    for i=1:size(resmd,1)
        resci = motion_resmd(i,c:c+corrwsize) ; 
        rescmts = mts(c:c+corrwsize) ; 
        xcorrs(count,i,:) = xcorr(rescmts,resci,20,'coeff') ; 
    end
    count = count + 1 ; 
end
figure,imagesc(squeeze(mean(xcorrs,1)),[-.6,.6]) ; axis xy ;

%}






mts = mean(ts,1) ; mts = mts(1:ntrs) ; 
clear xcorrs 
for i=1:size(resmd,1)
   xcorrs(i,:) = xcorr(mts(35:end-35),resmd(i,35:end-35),20,'coeff') ; 
end






%{
%mchans = log(squeeze(allchans)) ;
%basechans = mchans - repmat(mean(mchans(:,:,:),2),[1,size(mchans,2),1]) ; 
figure,for i=1:32 ; subplot(4,8,i) ; imagesc(squeeze(basechans(i,:,:))') ; axis xy ; end
figure,imagesc(squeeze(mean(basechans,1))') ; axis xy ; 
mbase = squeeze(mean(basechans,1))' ; 
[U,S,V] = svd(mbase) ;
f10 = (squeeze(basechans(:,:,10))) ; 
[U,S,V] = svd(f10) ; 
%}
%{

resbase = reshape(allchans,[64,size(allchans,2)*size(allchans,3)]) ; 
[weights,sphere] = runica(resbase,'maxsteps',128) ; 

acts = weights*sphere*resbase ; 
resacts = reshape(acts,size(allchans)) ; 

bads = zeros(1,64) ; 
for i=1:64 ; 
    subplot(5,13,i) ; 
    vec = (squeeze(mean(resacts(i,:,:),3))') ; 
    plot(vec) ; title(skewness(vec)) ; 
    if skewness(vec) > 4
        bads(i) = 1 ; 
    end
end
bads = find(bads==1) ;
for i=1:length(bads)
    subplot(5,7,i)
    vec = (squeeze(mean(resacts(bads(i),:,:),3))') ; 
    plot(vec) ; 
end

acts(bads,:) = 0 ; 
winv = pinv(weights*sphere) ; 
invacts = winv*acts ; 
resinvacts = reshape(invacts,size(allchans)) ;

%mchans = log(squeeze(resinvacts)) ;
%basechans = resinvacts - repmat(mean(resinvacts(:,:,:),2),[1,size(resinvacts,2),1]) ;
basechans = resinvacts ; 
for i=1:32 ; subplot(4,8,i) ; imagesc(log(abs(squeeze(allchans(i,:,:))'))) ; axis xy ; end
goodchans = [2,5,6,7,13] ;
mc = squeeze(mean(basechans(goodchans,:,:)))' ; 
mcprev = squeeze(mean(allchans(goodchans,:,:)))' ; 
mcbads = squeeze(mean(allchans(bads,:,:),1))' ;
plot(smooth(mean(mc(10:18,:)),100)) ; 

[s,f] = spectopo(ica2.icaact,0,250,'plot','off') ; 
figure,plot(mean(mc,2)) ; figure,plot(mean(mcbads,2)) ;
figure,imagesc(mc) ; axis xy ; figure,imagesc(mcbads) ; axis xy ; 
figure,plot(mean(mc(20:25,:),1)) ; hold on,plot(mean(mcprev(20:25,:),1),'r') ; 
figure,plot(mean(s(goodchans,:))) ; 
%}






