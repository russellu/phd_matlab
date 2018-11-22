clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
for sub=1%:length(subs)
cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls
mov = dir('bcgica*rest*set') ; 
EEG = pop_loadset(mov(1).name) ; 
EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
newz = EEG.data ; 
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
end

%subplot(3,3,i) ;
figure,imagesc(log(spec)',[0,8]) ; axis xy ;    
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

eeg2 = EEG ; eeg2.data = eegfiltfft(EEG.data,250,1,125) ;%+ eegfiltfft(EEG.data,250,11,15); 
icaw = load('icaw') ; icaw = icaw.icaw ; 
ica = eeg2 ; ica.icaweights = icaw{1} ; ica.icasphere = icaw{2} ; 
ica.icaact = ica.icaweights*ica.icasphere*EEG.data ; 
winv = pinv(ica.icaweights*ica.icasphere) ; for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),ica.chanlocs) ; title(i) ; end

M = 250 ; Wc = .1 ; k = (1:M-1); s = sin(Wc*k)./ k; c0 = [Wc,s]; A = toeplitz(c0);
[V,evals] = eig(A); % Only need the principal eigenvector
[emax,imax] = max(abs(diag(evals)));
w = V(:,end-6:end) ;clear g1 ; w = ones(1,length(M)) ; 
g1(:,1:size(w,2)) = w ; 
halfl = round(M/2) ; 
chans = 1:64 ; clear specs allchans 
for c=1:length(chans)
chan = ica.icaact(chans(c),:) ; icount = 1 ; 
for i=halfl+1:10:length(chan)-halfl
    datai = chan(i-halfl:i+halfl-1) ;
    windowed = repmat(datai',[1,size(w,2)]).*w ; 
    f = abs(fft(windowed,[],1)) ; 
    spec = f(1:halfl) ; 
    specs(icount,:) = spec ; 
    icount = icount + 1 ; 
end
allchans(c,:,:) = specs ; 
end

motion = (squeeze(mean(mean(allchans(:,:,1:25),1),3))') ; thresh = median(motion)*1.5 ; 
bads = find(motion>thresh) ; zmotion = zeros(size(motion)) ; zmotion(bads) = 1 ; 
zmotion = imdilate(zmotion,strel(ones(10,1))) ; 
bw = bwconncomp(zmotion==0) ; clusts = bw.PixelIdxList ; 
for i=1:length(clusts) ; if length(clusts{i})<10 ; zmotion(clusts{i}) = 1 ; end ; end
zmotion = imdilate(zmotion,strel(ones(4,1))) ;
mclusts = bwconncomp(zmotion) ; pclusts = mclusts.PixelIdxList ; 
mallchans = zeros(size(allchans)) ; 
for cc=1:size(allchans,1)
template = squeeze(allchans(cc,:,:)) ; 
for f=1:size(template,2)
for i=1:length(pclusts)
    
    indsi = pclusts{i} ; 
    if i==1 
        indsi(1:5) = [] ;
    elseif i==length(pclusts)
        indsi(end-5:end) = [] ;         
    end
    meanbefore = mean(template(indsi(1)-4:indsi(1),f)) ; 
    meanafter = mean(template(indsi(end):indsi(end)+4,f)) ; 
    template(indsi,f) = ones(1,length(indsi))*((meanbefore+meanafter)/2) ;
end
end
mallchans(cc,:,:) = template ; 
end


allchans = mallchans([13,14,16,19,21],:,:) ; 
mchans = (squeeze(mean(allchans,1))) ;
basechans = mchans - repmat(mean(mchans(:,:),1),[size(mchans,1),1]) ; 
md = medfilt2(basechans,[5,5]) ; 
md = eegfiltfft(md',250/10,0.02,10) ; md = md' ; 

evts = {EEG.urevent.type} ; 
lats = cell2mat({EEG.urevent.latency}) ; 
s1s = find(strcmp('S  1',evts)) ; 
r128s = find(strcmp('R128',evts)) ; 
s1lats = lats(s1s) ; r128lats = lats(r128s) ; 
windowlats = s1lats/10 ; windowr128s = r128lats/10 ; 
ntrs = length(r128s)-1 ; 

st = windowr128s(1) ; en = windowr128s(end) ; 
imagesc(md(st:en,:)',[-1,1]) ; axis xy ;% vline(windowlats,'k') ; 
newmd = md(st:en,:)' ; resmd = imresize(newmd,[halfl,ntrs]) ; 

cd(['c:/shared/newbadger_mri/',subs{sub}]) ; ls
bp = load_untouch_nii('bp_reg_topup_mc_retino_rest.nii.gz') ; bpimg = bp.img ; 

corr = load_untouch_nii('corrs.nii.gz') ; 
mask = medfilt3(corr.img > .3) ; inds = find(mask==1) ; [i1,i2,i3] = ind2sub(size(mask),inds) ; 
for i=1:length(inds)
    ts(i,:) = squeeze(bpimg(i1(i),i2(i),i3(i),:)) ; 
end
mts = mean(ts,1) ; mts = mts(1:ntrs) ; 




clear xcorrs 
wincr = 10 ; wlength = 80 ; wcount = 1 ; 
for ws=20:wincr:size(mts,2)-(wlength+35)

    for i=1:size(resmd,1)    
        fmri = mts(ws:ws+wlength) ; 
        eeg = smooth(resmd(i,ws:ws+wlength)) ;
       xcorrs(wcount,i,:) = xcorr(fmri,eeg,20,'coeff') ; 
    end
    wcount = wcount + 1 ; 
end


figure,imagesc(xcorrs,[-.8,.8]) ; axis xy ; 
save('xcorrs','xcorrs') ; 


%{
for s=1:length(subs)
   cd(['c:/shared/newbadger_mri/',subs{s}]) ; ls
    xcorrs = load('xcorrs') ; xcorrs = xcorrs.xcorrs ; 
    allxcorrs(s,:,:) = xcorrs ; 
end
%}





