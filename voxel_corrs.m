clear all ; close all 
comps = {[32,6,48],[9,18,20],[12,10,4],[24,9,16],[8,15,52],[11,21,33],[46,33,15],[25,16,21],[17,7,21]} ;

subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
subns = [1,2,3,4,6,7,8,9] ; 
for sub=3%:length(subns) ; 
cd(['c:/shared/badger_eeg2/',subs{subns(sub)}]) ; ls
mov = dir('bcgica*movie*set') ; 
EEG = pop_loadset(mov(1).name) ; 
EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
allica = load('highfreqs') ; allica = allica.highfreqs ; 
weights = allica{1} ; sphere = allica{2} ; EEG.icaact=weights*sphere*EEG.data ; 
ica = EEG ; ica.icaact = weights*sphere*EEG.data ; 
winv = pinv(weights*sphere) ; % figure,for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),EEG.chanlocs) ; title(i) ; end;
rawts = mean(EEG.data([61,29],:),1) - mean(EEG.data([62,30],:),1) ;

[s,f] = spectopo(ica.icaact(:,10000:end-10000),0,250,'plot','off') ;
meanspecs(sub,:) = (mean(s(comps{subns(sub)}(:),:),1)) ;

clear corrs
for i=1:64
    corrs(i) = corr2(smooth(abs(EEG.icaact(i,:)),250),smooth(abs(rawts),250)) ; 
end
[~,si] = sort(corrs,'descend') ; acts = ica.icaact(si(1:10),:) ; 
clear corrs ; 
corrdatas = zeros(size(EEG.data)) ; 
for i=1:size(EEG.data,1)
    corrs = corr(EEG.data(i,:)',EEG.data') ; 
    [~,si] = sort(corrs,'ascend') ; 
    corrdatas(i,:) = EEG.data(i,:) - EEG.data(si(1),:) ; 
end
sumdatas = sum(diff(corrdatas,1,2).^2,2) ; [sv,si] = sort(sumdatas,'descend') ; 
acts(11:20,:) = corrdatas(si(1:10),:) ; acts(1,:) = std(EEG.data,0,1) ; acts(2,:) = std(ica.icaact,0,1) ; 
acts(21:64,:) = corrdatas(11:54,:) ; 

freqs = 1:100 ; 
filts = zeros(size(ica.icaact,1),length(freqs),size(ica.icaact,2)) ; 
for j=1:length(freqs)
    filts(:,j,:) = eegfiltfft(ica.icaact,250,freqs(j)-1,freqs(j)+1) ; disp(j) ; 
end


lats = cell2mat({EEG.urevent.latency}) ; 
labs = {EEG.urevent.type} ; 
r128s = find(strcmp('R128',labs)) ; 
r128lats = round(lats(r128s)) ; 

cd(['c:/shared/newbadger_mri/',subs{subns(sub)}]) ;
fmri = load_untouch_nii('warp_movie.nii.gz') ; 



clear allres 
for i=1:64
    spec = squeeze(abs(filts(i,:,:))) ; disp(i) ; 
    tracts = spec(:,r128lats(1):r128lats(end)-20) ; 
    reshape_sz = length(r128s)-1 ; 
    restracts = imresize(tracts,[size(tracts,1),reshape_sz]) ;
    allres(i,:,:) = restracts ; 
end
motions = squeeze(allres(64,:,:))>0 ; tmotions = sum(motions,1) ;

chanres = squeeze(mean(allres,1)) ; 
mres = squeeze(mean(mean(allres(:,40:end,:),1),2)) ; 
motion = chanres ; 
binmotion = zeros(size(motion)) ; 
%allmotions(sub,:,:) = motion ; 
for i=1:size(motion,1)
    [z,mu,sigma] = zscore(motion(i,:)) ; 
    binmotion(i,:) = motion(i,:) > mu + sigma/2 ; 
end

bads = find(binmotion==1) ; goods = find(binmotion==0) ; 
[gx,gy] = meshgrid(1:size(binmotion,2),1:size(binmotion,1)) ; 

allbinmotions(sub,:,:) = binmotion ; 
newchans = allres ; 
for i=1:64 ; disp(i)  ;
    current_chan = squeeze(allres(i,:,:)) ; 
    badx = gx(bads) ; bady = gy(bads) ; goodx = gx(goods) ; goody = gy(goods) ; 
    vq = griddata(double(goodx),double(goody),double(current_chan(goods)),double(badx),double(bady),'cubic') ; 
    zres = current_chan ; zres(bads) = vq ; 
    newchans(i,:,:) = imfilter(zres,fspecial('gaussian',[3,9],3)) ; 
    rawchans(i,:,:) = zres ;
    rawmotionchans(i,:,:) = current_chan ; 
end

circshifts = 6:9 ; 
voxels = reshape(fmri.img,[numel(fmri.img(:,:,:,1)),size(fmri.img,4)]) ; 
mchans = squeeze(mean(newchans(comps{subns(sub)},:,:),1)) ; 
mraw = log(squeeze(mean(rawchans(comps{subns(sub)},:,:),1))) ; 
mrawmotion = log(squeeze(mean(rawmotionchans(comps{subns(sub)},:,:),1))) ; 
mraw(isnan(mraw)) = 0 ; mrawmotion(isnan(mrawmotion)) = 0 ; 
figure, subplot(1,2,1) ; imagesc((mrawmotion),[-4,1.5]) ; axis xy ;   xlabel('TR') ; ylabel('freq(hz)');subplot(1,2,2) ; imagesc((mraw),[-4,1.5]) ; axis xy ;

mchans(isnan(mchans)) = 0 ; 
voxels = voxels(:,1:size(mchans,2)+1) ; 

corrs = zeros(length(circshifts),100,size(voxels,1)) ; 
for i=1:length(circshifts) ; 
    shifty = circshift(mchans',circshifts(i)) ; 
    shifty = shifty' ;
    for j=1:100 ; disp(j) ; 
       corrs(i,j,:) = corr(shifty(j,20:end-20)',voxels(:,20:end-21)') ;  
    end
    
end
rescorrs = zeros(size(fmri.img,1),size(fmri.img,2),size(fmri.img,3),size(corrs,1),size(corrs,2)) ; 
for i=1:size(corrs,1)
    for j=1:size(corrs,2)
        rescorrs(:,:,:,i,j) = reshape(squeeze(corrs(i,j,:)),[size(fmri.img,1),size(fmri.img,2),size(fmri.img,3)]) ; 
    end
end

mrescorrs = squeeze(mean(rescorrs,4)) ; 
cat = load_untouch_nii('cat.nii.gz') ; 
%cat.img = mrescorrs ; save_untouch_nii(cat,'gamma1_corrs.nii.gz') ; 
end
