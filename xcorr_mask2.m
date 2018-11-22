clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
comps = {[14,21,20,27,19,26,16,12,4,10],[6,14,17,16,25,11,9,10],[6,10,16,28,27,11,13,15],[6,23,9,8,24,18],[12,18,15,17,36],[15,20,8,11,14,9],[6,21,15,10,22,14,12],[7,31,21,13,19,11,25],[8,14,18,28,19]} ;

for sub = 6:7%length(subs) 
cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls
mov = dir('bcgica*movie*set') ; 
EEG = pop_loadset(mov(1).name) ; 
rawts = mean(EEG.data([61,29],:),1) - mean(EEG.data([62,30],:),1) ;
EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
%EEG = pop_runica(EEG,'runica') ; tp(EEG) ; 
%mot = getmotion(EEG) ; bads = find(mot==0) ; goods = find(mot~=0) ; 
%EEG.icaact(64,:) = mot ; 
allica = load('allica') ; allica = allica.allica ; 
weights = allica{1} ; sphere = allica{2} ; EEG.icaact=weights*sphere*EEG.data ; tp(EEG) ; 


clear corrs
for i=1:64
    corrs(i) = corr2(smooth(abs(EEG.icaact(i,:)),250),smooth(abs(rawts),250)) ; 
end
clear specs ;
wincr = 25 ; 
[sv,si] = sort(corrs,'descend') ; 
for i=1:5 ; 
    specs(i,:,:) = compute_spectrogram(EEG.icaact(si(i),:),250,wincr) ; 
end
logspecs = log(specs) ; 
baselogs = logspecs - repmat(mean(logspecs,2),[1,size(logspecs,2),1]) ; 
%figure,imagesc(squeeze(mean(baselogs,1))') ; axis xy ; 
mbaselogs = squeeze(mean(baselogs,1)) ; zbase = zeros(size(mbaselogs)) ; 

[sv,si] = sort(mbaselogs(:),'ascend') ; 
zbaselogs = zeros(size(mbaselogs)) ; 
zbaselogs(si(1:round(length(si)/1.5))) = 1 ;

for i=1:125
    basei = mbaselogs(:,i).^2 ;
    badsi = find(basei>median(basei)*5) ; 
    zbase(badsi,i) = 1 ; 
end
bws = bwconncomp(zbase) ; szs = cellfun(@length,bws.PixelIdxList) ; 
bads = find(szs<25) ; badinds = bws.PixelIdxList(bads) ; for i=1:length(badinds) ; zbase(badinds{i}) = 0 ; end
%figure,imagesc(zbase') ; axis xy ; 
dil = (imdilate(zbase,strel(ones(5,3)))') ; dil = dil' ;% figure,imagesc(dil) ;
dil = (dil - zbaselogs) > 0 ; 

clear specs  ;
for i=1:size(EEG.icaact,1) ; disp(i) ; 
    specs(i,:,:) = compute_spectrogram(EEG.icaact(i,:),250,wincr) ;  
end
gridspecs = zeros(size(specs)) ; 
for i=1:64 ;  disp(i) ; 
[gx,gy] = meshgrid(1:size(dil,2),1:size(dil,1)) ;
goods = find(dil==0) ; goodx = gx(goods) ; goody = gy(goods) ; 
bads = find(dil==1) ; badx = gx(bads) ; bady = gy(bads) ; 
currspec = squeeze(specs(i,:,:)) ; 
goodvals = currspec(goods) ; 
vq = griddata(double(goodx),double(goody),double(goodvals),double(badx),double(bady),'linear') ; 
currspec(bads) = vq ; 
gridspecs(i,:,:) = currspec ; 

end


lats = cell2mat({EEG.urevent.latency}) ; 
labs = {EEG.urevent.type} ; 
r128s = find(strcmp('R128',labs)) ; 
r128lats = round(lats(r128s)/wincr) ; 

cd(['c:/shared/newbadger_mri/',subs{sub}])
corrs = load_untouch_nii('corrs.nii.gz') ; 
rest = load_untouch_nii('bp_reg_topup_mc_retino_movie.nii.gz') ; 
corrvox = find(medfilt3(corrs.img)>.12) ; [i1,i2,i3] = ind2sub(size(corrs.img),corrvox) ; ts = zeros(length(i1),size(rest.img,4)) ; 
for i=1:length(i1) ; ts(i,:) = squeeze(rest.img(i1(i),i2(i),i3(i),:)) ; end
mts = mean(ts(:,1:size(r128lats,2)-1)) ; 

clear allres ; 
for i=1:64
spec = squeeze(gridspecs(i,:,:)) ; disp(i) ; 
tracts = spec(r128lats(1):r128lats(end)-20,:) ; 
reshape_sz = length(r128s)-1 ; 
restracts = imresize(tracts,[reshape_sz,size(tracts,2)]) ;
%if i~=64
%restracts = eegfiltfft(restracts',1/0.693,0.01,1.5)' ; %restracts = imfilter(restracts,fspecial('gaussian',[5,5],1)) ;
%end
allres(i,:,:) = restracts ; 
end
motions = squeeze(allres(64,:,:))>0 ; tmotions = sum(motions,2) ;

clear allcorrs allbadinds
twincr = 100 ;
tstep = 20 ; 
for e=1:64 ; disp(e) ; 
    for j=1:125 ; tcount = 1 ; 
        for t=20:tstep:size(allres,2)-(twincr+20)
            c1 = mts(t:t+twincr) ; 
            s1 = squeeze(allres(e,t:t+twincr,j)) ; 
            newtmotions = tmotions(t:t+twincr) ; 
            xcoeff = masked_xcorr(c1,s1,newtmotions) ; 
            allcorrs(tcount,e,j,:) = xcoeff ; 
            tcount = tcount + 1 ; 
        end
    end
end
figure,for i=1:64 ; subplot(5,13,i) ; imagesc((squeeze(mean(allcorrs(:,i,:,:),1)))) ; title(i) ; axis xy ; end
ch = comps{sub}(1:2) ;
for i=1:30 ; subplot(5,6,i) ; imagesc(squeeze(mean(allcorrs(i,ch,:,:),2))) ; axis xy ; end
figure,
subplot(1,2,1) ; imagesc(squeeze(mean(mean(allcorrs(:,ch,:,:),1),2)),[-.2,.2]) ; vline(30,'k') ; axis xy
%subplot(2,2,2) ; plot(squeeze(mean(mean(mean(allcorrs(:,ch,12:30,:),1),2),3))) ; hline(0,'k') ; vline(30,'r') ; 
%subplot(2,2,4) ; plot((squeeze(mean(mean(mean(allcorrs(:,ch,:,30:45),1),2),3)))) ; hline(0,'k') ;
moviecorrs = allcorrs ; save('moviecorrs','moviecorrs','-v7.3') ; 
end
