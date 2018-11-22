clear all ; close all 
occ = {[65,13],[44,41,28],[60,48],[85,82,67,6],[81,54,36],[77,76,74],[53,49],[66,53,48,38],[77,75,42]} ; 
lingual = {[80,16],[8],[18],[50,5],[60],[53],[31],[6],[16]} ; 
lateral = {[17],[23],[41],[21],[19],[51],[33],[24],[22]} ; 
dmn = {[15,10],[7,2],[12,9,8],[8,1],[1],[33],[43,38],[12,7],[5,3,2]} ; 
motor = {[71,22],[20],[21],[32,11],[53,51,10],[31],[40,5],[88,71,49,8],[51,21,15]} ; 
muscle = {[85],[],[65],[63,61],[85],[87],[40],[60],[81]} ; 
csf = {[5],[5],[3],[31],[48],[6],[29],[33],[29]} ; 
motion = {[33,29],[12],[57],[76,72,30],[59,27],  [34],[63],[81,45],[66,44,43]} ; 
eyes = {[76],[],[70],[80],[26],[55],[55],[51],[41]} ; 
misc_cortex = {[56,54,53,39,35,32,24,9],[34,25,23,15],[30,28,26,31,19,3],[65,55,49,25,22,20,16,14,13,9,4],...
[45,32,28,20,17,16,15,11,6,5,4],[69,68,60,59,29,25,20,18,17,15,13,12],[33,30,12],[72,61,52,47,42,35,34,80,23,22,14,11],[64,33,32,27,24,13,11,7]} ; 
white = {[51],[10],[27],[19],[8],[36],[20],[17],[30]} ; 
parietal = {[49],[7],[36],[34,12,7],[9],[71,9],[38],[20],[12,6]} ; 
cbllm = {[31],[18],[58],[75],[77],[26],[26],[29],[36]} ; 
back = {[1],[21],[36],[71],[69],[16],[58],[19,15],[18]} ;  


% EEG comps 
comps = {[14,21,20,27,19,26,16,12,4,10],[6,14,17,16,25,11,9,10],[6,10,16,28,27,11,13,15],[6,23,9,8,24,18],[12,18,15,17,36],[15,20,8,11,14,9],[6,21,15,10,22,14,12],[7,31,21,13,19,11,25],[8,14,18,28,19]} ;

subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
subns = [1,2,3,4,6,7,8,9] ; 
for sub=1%:length(subns) ; 
cd(['c:/shared/badger_eeg2/',subs{subns(sub)}]) ; ls
mov = dir('bcgica*movie*set') ; 
EEG = pop_loadset(mov(1).name) ; 
EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
%ica = pop_runica(EEG,'runica') ;
allica = load('allica') ; allica = allica.allica ; 
weights = allica{1} ; sphere = allica{2} ; EEG.icaact=weights*sphere*EEG.data ; 
ica = EEG ; ica.icaact = weights*sphere*EEG.data ; 
winv = pinv(weights*sphere) ; figure,for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),EEG.chanlocs) ; title(i) ; end;
rawts = mean(EEG.data([61,29],:),1) - mean(EEG.data([62,30],:),1) ;
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
wincr = 50 ; 

clear specs 
for i=1:64 ; disp(i) ; 
   specs(i,:,:) = compute_spectrogram(acts(i,:),250,wincr) ;  
end
logspecs = log(specs) ; 
baselogs = logspecs - repmat(mean(logspecs,2),[1,size(logspecs,2),1]) ; 
figure,imagesc((squeeze(mean(baselogs(1:20,:,:),1)))') ; axis xy ; 
mbaselogs = squeeze(mean(baselogs(1:20,:,:),1)) ; zbase = zeros(size(mbaselogs)) ; 
[s,v,d] = svd(mbaselogs') ; 
repabs = repmat(abs(d(:,1)),[1,size(mbaselogs,2)]) ; 
figure,imagesc((repabs.*mbaselogs)') ; 
mat = abs(repabs.*mbaselogs) ; binmat=  zeros(10,size(mat)) ; 
for mthresh=1:10
for i=1:size(mat,2)
    medi = median(mat(:,i)) ; 
    binmat(mthresh,:,i) = mat(:,i) > (medi*mthresh) ; 
end
end
clear icaspecs 
for i=1:size(ica.icaact,1) ; disp(i) ; 
    icaspecs(i,:,:) = compute_spectrogram(ica.icaact(i,:),250,wincr) ; 
end


sumbin = squeeze(sum(binmat,1)) ; mlogs = squeeze(mean(logspecs,1)) ; 
thresh5 = sumbin>2 ; [gx,gy] = meshgrid(1:size(sumbin,1),1:size(sumbin,2)) ; figure,imagesc(thresh5') ; axis xy ;

clear allcorrected
for i=1:64 ; disp(i) ; 
    current_specs = squeeze(icaspecs(i,:,:)) ; 
    badinds = find(thresh5==1) ; goodinds = find(thresh5==0) ; 
    goodx = gx(goodinds) ; goody = gy(goodinds) ; goodv = current_specs(goodinds) ; 
    badx = gx(badinds) ; bady = gy(badinds) ; 
    vq = griddata(double(goodx),double(goody),double(goodv),double(badx),double(bady)) ; 
    current_specs(badinds) = vq ; 
    allcorrected(i,:,:) = current_specs ; 
end

lats = cell2mat({EEG.urevent.latency}) ; 
labs = {EEG.urevent.type} ; 
r128s = find(strcmp('R128',labs)) ; 
r128lats = round(lats(r128s)/wincr) ; 
%{
cd(['c:/shared/newbadger_mri/',subs{sub}])
corrs = load_untouch_nii('corrs.nii.gz') ; 
rest = load_untouch_nii('bp_reg_topup_mc_retino_movie.nii.gz') ; 
corrvox = find(medfilt3(corrs.img)>.12) ; [i1,i2,i3] = ind2sub(size(corrs.img),corrvox) ; ts = zeros(length(i1),size(rest.img,4)) ; 
for i=1:length(i1) ; ts(i,:) = squeeze(rest.img(i1(i),i2(i),i3(i),:)) ; end
mts = mean(ts(:,1:size(r128lats,2)-1)) ; 
%}
cd(['c:/shared/newbadger_mri/',subs{subns(sub)},'/melodic']) ;
mix = load('melodic_mix') ; 
filtmix = eegfiltfft(mix',1/0.693,0.02,1.5)' ; 
restmix = filtmix(end-449:end,:) ; 
moviemix = filtmix(end-(449+735):end-450,:) ; 
retinomix1 = filtmix(1:735,:) ; 
retinomix2 = filtmix(736:736+734,:) ; 
gammamix1 = filtmix(1471:1471+734,:) ; 
gammamix2 = filtmix(2206:2206+734,:) ; 

ots = mean(moviemix(:,occ{subns(sub)}),2) ; 
mts = ots(1:end-1) ; 

%chans = [6,10,11] ; 
%allcorrected(1,:,:) = squeeze(mean(allcorrected(chans,:,:),1)) ; 
clear allres 
for i=1:64
spec = squeeze(allcorrected(i,:,:)) ; disp(i) ; 
tracts = spec(r128lats(1):r128lats(end)-20,:) ; 
reshape_sz = length(r128s)-1 ; 
restracts = imresize(tracts,[reshape_sz,size(tracts,2)]) ;
if i~=64
%restracts = eegfiltfft(restracts',1/0.693,0.01,1.5)' ; %restracts = imfilter(restracts,fspecial('gaussian',[5,5],1)) ;
end
allres(i,:,:) = restracts ; 
end
motions = squeeze(allres(64,:,:))>0 ; tmotions = sum(motions,2) ;
clear allcorrs allbadinds
for e=1:64 ; disp(e) ; 
    for j=1:125 ; tcount = 1 ; 
        c1 = mts(20:end-20) ; 
        s1 = squeeze(allres(e,20:end-20,j)) ; 
        newtmotions = tmotions(20:end-20) ; 
        xcoeff = masked_xcorr(c1,s1,newtmotions) ; 
        allcorrs(e,j,:) = xcoeff ; 
    end
end
allsubcorrs(sub,:,:,:) = allcorrs ; 
%figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(allcorrs(i,:,:))) ; title(i) ; axis xy ; end
%chans = [8,15,18,24] ; 
%figure,imagesc((squeeze(mean(allcorrs(chans,:,:),1))),[-.20,.20]) ; axis xy
end

for i=1:8 ; figure ; for j=1:64 ; subplot(5,13,j) ; imagesc(squeeze(allsubcorrs(i,j,:,:))) ; axis xy ; end ; end



for i=1:8 ; 
    subplot(3,3,i) ;
    imagesc(squeeze(mean(allsubcorrs(i,comps{subns(i)}(1:3),:,:),2)),[-.2,.2]) ; axis xy ; 
    means(i,:,:) = squeeze(mean(allsubcorrs(i,comps{subns(i)}(1:2),:,:),2)) ; 
end
figure,
subplot(1,2,1) ; imagesc(squeeze(mean(means,1))) ; axis xy ;
subplot(2,2,2) ; shadedErrorBar([],squeeze(mean(mean(means(:,10:25,:),2),1)),squeeze(std(mean(means(:,10:25,:),2),0,1))/sqrt(8)) ; hline(0,'k') ; vline(30,'r') ; hold on ; 
shadedErrorBar([],squeeze(mean(mean(means(:,40:90,:),2),1)),squeeze(std(mean(means(:,40:90,:),2),0,1))/sqrt(8),{'r'})
subplot(2,2,4) ; shadedErrorBar([],squeeze(mean(mean(means(:,:,35:40),3),1)),squeeze(std(mean(means(:,:,35:40),3),0,1))/sqrt(8)) ; hline(0,'k') ; 



for i=1:125
    for j=1:60
        [h,p,ci,stats] = ttest(squeeze(means(:,i,j))) ; 
        ts(i,j) = stats.tstat ; 
    end
end

figure,imagesc(ts,[-6,6]) ; axis xy ; colorbar ; 





