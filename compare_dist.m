clear all  ;close all
subs = {'alex','genevieve','dina','jeremie','russell','tegan','valerie'} ; 
for s=1:length(subs)
cd(['C:\shared\simdenoise\',subs{s}]) ; 
bdersp = load('bdersp') ; bdersp = bdersp.bdersp ; 
EEG = pop_loadbv('.','retino_gamma_02.vhdr') ; 
EEG = pop_chanedit(EEG,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ; 

cd(['C:\shared\all_white_normals\a2_good\sub1_',subs{s}]) ;
corrs = load_untouch_nii('meancorrs_padt1.nii.gz') ; 
allmasks(:,:,:,s) = corrs.img;  
locs = load_untouch_nii('t1_coords.nii.gz') ; 
length(unique(locs.img(:)))
clear cx cy cz
for i=1:65
   [cx(i),cy(i),cz(i)] = centmass3(locs.img==i) ;  
    
end
[rx,ry,rz] = centmass3(corrs.img>.5) ; 
order = load('c:/shared/all_white_normals/a2_good/elecorder.mat') ; order = order.elecorder ; 
labs = {EEG.chanlocs.labels} ; clear allinds 
for i=1:length(labs)
    ind = find(strcmpi(labs{i},order)) ; 
    if i==32
    allinds(i) = 0 ; 
    else allinds(i) = ind ; 
    end
end ; allinds(32) = 1 ; 

melecs = squeeze(mean(mean(bdersp(1,:,:,:,40:150),3),5)) ; 
allelecs(s,:,:) = melecs ; 
mdx = cx(allinds) ; mdy = cy(allinds) ; mdz = cz(allinds) ; 

allcoords(s,1,:) = mdx ; allcoords(s,2,:) = mdy ; allcoords(s,3,:) = mdz ; 

diffs = [mdx-rx;mdy-ry;mdz-rz] ;
sqrdiff = sqrt(sum(diffs.^2,1)) ; 

sqrdiff(32) = [] ; 
melecs(32,:) = [] ; 
alldiffs(s,:) = sqrdiff ; %allelecs(s,:,:) = melecs ; 
c = corr(sqrdiff',melecs) ;
cs(s,:) = c ; 
end

% make some plots for the small group 
%plot(mean(alldiffs([1,2,5,6,7],:),1),squeeze(mean(allelecs([1,2,5,6,7],:,29),1)),'o'); hline(0,'k') ; 

% try a different distance metric, instead of center of mass only: 
% maybe mean distance to 100 closest voxels or something

for i=1:7
    maski = allmasks(:,:,:,i) > 0.25 ; maski(:,1:180,:) = 0 ; 
   % subplot(2,4,i) ; imagesc(sum(maski,3)>0) ; 
    for j=1:64
        maskinds = find(maski==1) ; 
        [maskx,masky,maskz] = ind2sub(size(maski),maskinds) ; 
        coordsij = squeeze(allcoords(i,:,j)) ; 
        diffs = [maskx-coordsij(1),masky-coordsij(2),maskz-coordsij(3)] ; 
        sqrdiffs = sqrt(sum(diffs.^2,2)) ; 
        [sv,si] = sort(sqrdiffs,'ascend') ; 
        alldists(i,j) = mean(sv) ; 
    end
end
newelecs = allelecs ; newdists = alldists ; 
newelecs(:,32,:) = [] ; newdists(:,32,:) = [] ; 
clear corrs ; 
for i=1:7
    corrs(i,:) = corr(newdists(i,:)',squeeze(newelecs(i,:,:))) ; 
end

shadedErrorBar([],mean(corrs,1),std(corrs,1)/sqrt(7)) ; hline(0,'k') ; hold on ; 
for i=1:60 ; [h,p,ci,stats] = ttest(corrs(:,i)) ; allts(i) = stats.tstat ; allps(i) = p ; end ; 
for i=1:60 ; if allps(i) < 0.05/60 ; text(i,.55,'*') ; end ; end ; set(gca,'XTick',0:5:60,'XTickLabel',0:10:120) ; xlabel('frequency(hz)') ; ylabel('corr(rho)') ; 
title('mean distance vs power correlations (n=7), *p<0.01,uncorrected') ; 
figure,plot(mean(newdists,1),squeeze(mean(mean(newelecs(:,:,25:30),1),3)),'o') ; xlabel('distance(mm)') ; ylabel('gamma modulation (db)') ; 

for i=1:63 ; scorrs(i,:) = corr(squeeze(newelecs(:,i,:)),newdists(:,i)) ; end
newcorrs = zeros(64,60) ; newcorrs(1:31,:) = scorrs(1:31,:) ; newcorrs(33:end,:) = scorrs(32:end,:) ; 
topoplot(mean(newcorrs(:,25:30),2),EEG.chanlocs,'maplimits',[-.8,.8]) ; 

plot(mean(alldists(:,[9,20,10]),2),squeeze(mean(mean(allelecs(:,[9,20,10],28:31),2),3)),'o') ; lsline ; xlabel('distance(mm)') ; ylabel('power(db)') ; 

subplot(1,2,1) ; topoplot(squeeze(mean(mean(allelecs(:,:,5:10),1),3)),EEG.chanlocs,'electrodes','off','maplimits',[-1.5,1.5]);
subplot(1,2,2) ; topoplot(squeeze(mean(mean(allelecs(:,:,20:30),1),3)),EEG.chanlocs,'electrodes','off','maplimits',[-1.5,1.5]);


