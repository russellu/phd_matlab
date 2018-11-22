clear all ; close all ; 
subs = {'charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','lisa','marc','marie',...
    'mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ;

cd c:/shared/all_white_normals/a2_good/ ; 
elabs = load('elabs.mat') ; elecorder = load('elecorder.mat') ; 
elabs = elabs.elabs ; elecorder = elecorder.elecorder ;

clear allersp ; 
for s=1:length(subs)
    cd(['c:/shared/allres/',subs{s}]) ; ls 
    ersp = load('ersp.mat') ; 
    stersp = load('stersp') ; stersp = stersp.ersp ; 
    allstersp(s,:,:,:,:) = squeeze(stersp(:,1:3,:,:)) ;
    ersp = ersp.ersp ; allersp(s,:,:) = squeeze(mean(ersp(1:2,:,:),1)) ; 
  
    cd(['c:/shared/all_white_normals/a1_good/sub_',subs{s}]) ; 
    corrs = load_untouch_nii('noeyecorrs.nii.gz') ; 
    elecbrain = load_untouch_nii('elecbrain.nii.gz') ; 

    leftnorms = dir('*whitesub*lh*gz') ; rightnorms = dir('*whitesub*rh*gz') ; 
    leftnorms = load_untouch_nii(leftnorms.name) ; rightnorms = load_untouch_nii(rightnorms.name) ; 
    leftbin = sum(leftnorms.img>0,4)>0 ; rightbin = sum(rightnorms.img>0,4)>0 ; 
    
    coords = load('coords') ; coords = coords.coords ; 
    bincorrs = corrs.img > .35 ; 
    bincorrs(:,1:end-120,:) = 0 ; 
    leftcorrs = leftbin.*bincorrs ; rightcorrs = rightbin.*bincorrs ; 
    
    figure,imagesc(squeeze(sum(bincorrs,3))) ; 
    [cx,cy,cz] = centmass3(bincorrs) ; 
    [rcx,rcy,rcz] = centmass3(rightcorrs) ; 
    [lcx,lcy,lcz] = centmass3(leftcorrs) ; 
    
    dcoords = [coords(1,:)-cx;coords(2,:)-cy;coords(3,:)-cz] ; 
    dcoords = sqrt(sum(dcoords.^2,1)) ;
    alldcoords(s,:) = dcoords ; 
    
    rdcoords = [coords(1,:)-rcx;coords(2,:)-rcy;coords(3,:)-rcz] ; 
    rdcoords = sqrt(sum(rdcoords.^2,1)) ;
    allrdcoords(s,:) = rdcoords ; 
    
    ldcoords = [coords(1,:)-lcx;coords(2,:)-lcy;coords(3,:)-lcz] ; 
    ldcoords = sqrt(sum(ldcoords.^2,1)) ;
    allldcoords(s,:) = ldcoords ; 
    
end

cd('c:/shared/allres/alex') ; 
EEG = pop_loadset('resamp_vis01.set') ; 
labs = {EEG.chanlocs.labels} ;
for i=1:length(labs)
    indi = find(strcmpi(labs{i},elecorder)) ; 
    if ~isempty(indi)
        allinds(i) = indi ; 
    end
end
allinds(28) = allinds(29) ; allinds(32) = allinds(31) ; 


mersp = squeeze(mean(mean(mean(allstersp(:,:,:,:,40:180),2),3),5)) ;
subplot(1,3,3) ; 
[corrmat,pmat] = corr(mersp,allrdcoords) ; 
corrvals = mean(corrmat(30:45,allinds),1) ; 
topoplot(corrvals,EEG.chanlocs) ; colorbar

subplot(1,3,2) ; 
corrmat = corr(mersp,alldcoords) ; 
corrvals = mean(corrmat(30:45,allinds),1) ; 
topoplot(corrvals,EEG.chanlocs) ; colorbar

subplot(1,3,1) ; 
corrmat = corr(mersp,allldcoords) ; 
corrvals = mean(corrmat(30:45,allinds),1) ; 
topoplot(corrvals,EEG.chanlocs) ; colorbar



plot(squeeze(allrdcoords(:,55)),squeeze(mersp(:,45)),'o') ; lsline ; 
title(['rho=',num2str(corrmat(45,55)),' p=',num2str(pmat(45,55))]) ; 
xlabel('distance electrode to cortex (mm)') ; ylabel('high gamma power (60-90Hz) modulation (dB)') ; 

plot(corrmat(:,55),'LineWidth',2) ; hold on ; plot(corrmat(:,54),'r','LineWidth',2) ;plot(corrmat(:,59),'m','LineWidth',2) ;
hline(0,'k') ; legend({elecorder{55},elecorder{54},elecorder{59}}) ; 
set(gca,'XTick',1:5:60,'XTickLabel',1:10:120) ; xlabel('frequency(hz)') ; ylabel('correlation (spearmans rho)') ;


nsubs = {'alex','alexandra3','audrey','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','leila','lisa','marc',...
    'marie','mathieu','maxime','menglu','mingham','olga','patricia','po','russell','suhan2','sunachakan','tah','tegan2','vincent'} ; 

for i=1:length(subs) ; subinds(i) = find(strcmpi(subs{i},nsubs)) ; end

fullcomps = {[11,12,6,13,30,38],[13,10,4,25,32,17,38,6,9],[2,5,7,8,9,15,16,22,26],[2,6,8,17,11],[5,7,6,8],[8,6,11,19,34,4],...
    [8,6,5,12],[10,14,16,22,28],[10,4,15,11],[9,5,11,12],[19,4,34],[5,16,3,26,17],[4,6,15,11,17,25],[6,13,22,29,54],[15,24,13],...
    [5,6,4,10,14],[8,1,5,13,14],[5,6,13,21],[1,13,10,18,22,23],[5,7,11,19],[4,10,5,6,8,9,29],[9,3,7,5,14],[7,6,3,17,25],[7,3,6,15,19],...
    [13,9,16],[20,18,17,3,5],[3,7,11,13,21],[8,9,5,18,44],[9,6,11,24],[9,6,8,13,15],[3,5,9,16]} ;
nfullcomps = fullcomps(subinds) ; 

for s=1:length(subs) 
   cd(['c:/shared/allres/',subs{s}]) ;  
   fullspec = load('fullspec.mat')  ;
   weights = fullspec.fullspec{1} ; sphere = fullspec.fullspec{2} ; 
   winv = pinv(weights*sphere) ; 
   rcomps(s,:) = winv(:,nfullcomps{s}(1)) ; 
   lcomps(s,:) = winv(:,nfullcomps{s}(2)) ; 
   mcomps(s,:) = winv(:,nfullcomps{s}(3)) ; 
   cd(['c:/shared/all_white_normals/a1_good/sub_',subs{s}]) ; ls 
   t1 = load_untouch_nii('t1.nii.gz') ; 
   corrs = load_untouch_nii('noeyecorrs.nii.gz') ; corrs.img(:,1:end-120,:) = 0 ; 
   bcorrs = corrs.img > .3 ; [cx,cy,cz] = centmass3(bcorrs) ; corrs.img = corrs.img.*double(corrs.img>.1) ; 
   tslices{s} = t1.img(:,:,cz) ; fslices{s} = mean(corrs.img(:,:,cz-5:cz+5),3) ; 
end
subplot(1,3,1) ; topoplot(mean(lcomps,1),EEG.chanlocs) ;
subplot(1,3,2) ; topoplot(mean(mcomps,1),EEG.chanlocs) ;
subplot(1,3,3) ; topoplot(mean(rcomps,1),EEG.chanlocs) ;

for i=1:length(tslices)
    subplottight(3,8,i) ; 
   plotoverlayIntensity2D(tslices{i},mat2gray(fslices{i}),fslices{i},270) ; title(['subject #',num2str(i)]) ; 
    
end

freqs = 1:2:120 ; times = -.5:.015:2.5 ; 
for i=1:22 ; subplot(4,6,i) ; 
    imagesc(times,freqs,squeeze(mean(mean(allstersp(i,:,:,:,:),2),3)),[-2,2]) ; 
    axis xy ; title(['subject #',num2str(i)]) ;
end
suptitle('inter-subject EEG variability (n=22)') ; 


for i=1:6 ; subplot(1,6,i)
    imagesc(times,freqs,squeeze(mean(mean(allstersp(:,i,:,:,:),1),3)),[-1.5,1.5]) ; axis xy ; 
    if i==1; xlabel('time(s)') ; ylabel('frequency(hz)') ; end
end




bar(squeeze(mean(mean(mean(mean(allstersp(:,:,:,5:12,40:180),2),3),4),5))); hold on ;
bar(squeeze(mean(mean(mean(mean(allstersp(:,:,:,20:45,40:180),2),3),4),5)),'r');
mfreqs = squeeze(mean(mean(mean(allstersp(:,:,:,:,40:180),2),3),5)) ; 

subplot(1,2,1) ; 
imagesc(freqs,1:22,mfreqs,[-2,2]) ; ylabel('subject #') ; xlabel('frequency(hz)') ; 
subplot(1,2,2) ; 
imagesc(freqs,freqs,c,[-.4,1]) ; xlabel('frequency(hz)') ; ylabel('frequency(hz)') ; 

for i=1:22 ; subplot(2,11,i) ; imagesc(times,freqs,squeeze(mean(mean(allstersp(i,:,:,:,:),2),3)),[-2,2]) ; axis xy ; title(['subject #',num2str(i)]) ; end

subplot(2,2,1) ; imagesc([-2,2]) ; colorbar ; subplot(2,2,2) ; imagesc([-.4,1]) ; colorbar ; subplot(2,2,3) ; imagesc([-1.5,1.5]) ; colorbar ;

