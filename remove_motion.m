clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
for sub=6%:length(subs)
cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls
mov = dir('bcgica*rest*set') ; 
EEG = pop_loadset(mov(1).name) ; 
EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;

data = EEG.data(:,32700:33000) ;
[u,s,v] = svd(data) ; 
u(:,1:4) = 0 ; 
sdata = u*s*v' ; 

%subplot(1,2,1) ; plot(sdata') ; ylim([-500,500]) ; subplot(1,2,2) ; plot(data') ; ylim([-500,500]) ;

%for i=1:64; subplot(5,13,i) ; topoplot(double(u(:,i)),EEG.chanlocs) ; title(i) ; end

wsize = 1000 ; icount = 1 ;  
newdata = zeros(size(EEG.data)) ; 
for i=1:wsize:size(EEG.data,2)-wsize
    data = EEG.data(:,i:i+wsize) ;
    [u,s,v] = svd(data) ; 
    topos(icount,:) = u(:,4) ; 
    
    icount = icount + 1 ;
    u(:,1:3) = 0 ; 
    sdata = u*s*v' ; 
    newdata(:,i:i+wsize) = sdata ; 
end
%for i=1:100 ; subplot(10,10,i) ; topoplot(double(topos(i,:)),EEG.chanlocs) ; end
[weights,sphere] = runica(newdata,'maxsteps',128) ; 
winv = pinv(weights*sphere) ; 
acts = weights*sphere*newdata ; 
figure,for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),EEG.chanlocs) ; title(i) ; end
[s,f] = spectopo(acts,0,250,'plot','off') ; 
end



M = 250 ; Wc = .1 ; k = (1:M-1); s = sin(Wc*k)./ k; c0 = [Wc,s]; A = toeplitz(c0);
[V,evals] = eig(A); % Only need the principal eigenvector
[emax,imax] = max(abs(diag(evals)));
%w = V(:,end-2:end) ;clear g1 
w = ones(1,length(M)) ; 
g1(:,1:size(w,2)) = w ; 
halfl = round(M/2) ; wincr = 25 ; 
chans = [1:64] ; clear specs allchans 
for c=1:length(chans) ; disp(['chan=',num2str(c)]) ; 
chan = acts(chans(c),:) ; icount = 1 ; 
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

tracts = basechans(:,r128lats(1):r128lats(end),:) ; 
reshape_sz = length(r128s)-1 ; restracts = zeros(size(tracts,1),reshape_sz,size(tracts,3)) ; 
for i=1:size(tracts,1) ; restracts(i,:,:) = imresize(medfilt2(squeeze(((tracts(i,:,:)))),[5,5]),[reshape_sz,size(tracts,3)]) ; end
filteeg = zeros(size(restracts)) ; 
for i=1:64 ; rts = eegfiltfft(squeeze(restracts(i,:,:))',1/0.693,0.02,1.5) ; filteeg(i,:,:) = rts' ; end


meanres = squeeze(mean(restracts(chans,:,:))) ; 
wincr = 10 ; wsize = 75 ; tcount = 1 ; clear xcorrs ; 
for t=35:wincr:size(mts,2)-(wsize+35)
    for i=1:64 ; 
        for j=1:125
            fmri = mts(t:t+wsize) ;
            eeg = squeeze(restracts(i,t:t+wsize,j)) ; 
            %eeg = meanres(t:t+wsize,j) ; 
            xcorrs(tcount,i,j,:) = xcorr(fmri,eeg,25,'coeff') ; 
        end
    end 
    tcount = tcount + 1 ; 
end

%tp(EEG) ; 
for d=1:31 ; 
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(xcorrs(d,i,:,:))) ; title(i) ; axis xy ; end
end


chans = [11] ;
imagesc(squeeze(mean(mean(xcorrs(:,chans,:,:),1),2))) ;  axis xy ; 


for i=1:31 ; subplot(5,7,i) ; imagesc(squeeze(mean(xcorrs(i,chans,:,:)))) ; axis xy ; title(i)  ;end ; 

gooddyns = [1:31] ; 
imagesc(squeeze(mean(mean(xcorrs(gooddyns,chans,:,:),1),2)),[-.3,.3]) ;  axis xy ; vline(25,'k') ; 







