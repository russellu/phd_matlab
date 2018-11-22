clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 
subgenders = [1,1,1,1,1,0,0,0,1,1,0,0,0,1,0,1,1,1,0,1,1,0,0,1]; 

peaksubs = [2,3,4,5,6,7,8,10,11,12,13,15,17,18,20,21,22,23] ; 
allpeaks = [64,76,74,64,66,84,72,66,74,74,70,66,64,70,64,88,70,56,84,60,66,66,54,84] ; 
stimtitles = {'unperturbed','5% contrast','33% contrast','plaid','10% random','60% random'} ; 

elabs = load('C:\shared\all_white_normals\a2_good\elabs') ; elabs = elabs.elabs ; 
elecorder = load('C:\shared\all_white_normals\a2_good\elecorder') ; elecorder = elecorder.elecorder ; 
baseeg = pop_loadset('c:/shared/allres/alex/cleanfilt.set') ; 
baselabs = {baseeg.chanlocs.labels} ;
for i=1:length(baselabs) 
    if ~isempty(find(strcmpi(baselabs{i},elecorder)))
        orderinds(i) = find(strcmpi(baselabs{i},elecorder)) ;
    end
end 
times = load('c:/shared/allres/alex/times') ; times = times.times ; 
freqs = load('c:/shared/allres/alex/freqs') ; freqs = freqs.freqs ; 
for sb=1:length(subs); disp(sb) ; 
    cd(['c:/shared/allfmris/sub_',subs{sb}]) ;  
    mstims = load('mstims') ; mstims = mstims.mstims; allmstims(sb,:,:) = mstims; 

    locs = load('locs') ; locs = locs.locs ; 
    meancorrs = load_untouch_nii('trim_cleancorrs_fs.nii.gz') ; t1 = load_untouch_nii('fs_t1.nii.gz') ; 
    
    gm = load_untouch_nii('gm.nii.gz'); 
    sumgms(sb) = sum(gm.img(:)); 
    
    threshs = 0.1:0.05:0.9 ; 
    for thr=1:length(threshs)
    threshcorrs = meancorrs.img>threshs(thr) ; sumcorrs(sb,thr) = sum(threshcorrs(:)) ; 
    end
    [sv,si] = sort(meancorrs.img(:),'descend') ; 
    bininds = find((meancorrs.img(:)>.6)); 
    zcorrs = zeros(size(meancorrs.img)) ; zcorrs(si(1:8000)) = 1 ; 
    [sx,sy,sz] = ind2sub(size(meancorrs.img),si(1:8000)) ; 
    [cx,cy,cz] = centmass3(zcorrs) ; 
    sqrdiffs = sqrt((locs(1,:)-cx).^2 + (locs(2,:)-cy).^2 + (locs(3,:)-cz).^2) ; 
    for i=1:size(locs,2)
        for j=1:length(sx)
            mindiffs(i,j) = sqrt((locs(1,i)-sx(j)).^2 + (locs(2,i)-sy(j)).^2 + (locs(3,i)-sz(j)).^2) ;
        end
    end
    mindiffs = min(mindiffs,[],2) ;
    distmat = zeros(1,64) ; 
    for i=1:length(orderinds) ; if orderinds(i)~=0 ; distmat(i) = sqrdiffs(orderinds(i)) ; end ; end
    for i=1:length(orderinds) ; if orderinds(i)~=0 ; mindistmat(i) = mindiffs(orderinds(i)) ; end ; end
    bads = find(distmat==0) ; goods = find(distmat~=0) ; 
    goods(find(goods==17)) = [] ; goods(find(goods==22)) = [] ; goods(find(goods==41)) = [] ; goods(find(goods==46)) = [] ; 
    alldists(sb,:) = distmat ; 
    cd(['c:/shared/allres/',subs{sb}]) ; 
    alldiffs = load('alldiffs.mat') ; alldiffs = alldiffs.alldiffs ; alleeg(sb,:,:,:) = alldiffs ; 
    
    %inot
    cd(['c:/shared/gamma_t1s/sub_',subs{sb},'/mri']) ;  
    leftnorms = load_untouch_nii('surf.nii.gz') ; 
    rightnorms = load_untouch_nii('surfr.nii.gz') ; 
    bothnorms = leftnorms.img + rightnorms.img ; 
    normls = 100:500:15000 ; 
    for nl=1:length(normls)
    [nx,ny,nz] = ind2sub(size(zcorrs),si(1:normls(nl))) ; 
    [bx,by,bz] = ind2sub(size(zcorrs),bininds); 
    allnorms = zeros(length(nx),3) ; 
    allbnorms = zeros(length(bx),3); 
    for i=1:length(nx)
        allnorms(i,:) = squeeze(bothnorms(nx(i),ny(i),nz(i),:)) ; 
    end
    for i=1:length(bx)
        allbnorms(i,:) = squeeze(bothnorms(bx(i),by(i),bz(i),:)) ; 
    end
    binots(sb) = 1 - sum(abs(sum(allbnorms,1)))./length(bx) ; 
    inots(sb,nl) = 1 - sum(abs(sum(allnorms,1)))./length(nx) ; 
    end
    cd(['c:/shared/allres/',subs{sb}]) ;  
    mersp = load('mersp') ; mersp = mersp.mersp ; allmersp(sb,:,:) = mersp ; 
    amersp = load('amersp') ; amersp = amersp.amersp ; allamersp(sb,:,:,:) = amersp ; 
end
bades = [41,17,28,32,22,46] ; goods = zeros(1,64) ; goods(bades)=1 ; goods = find(goods==0) ; 
for i=1:size(alleeg,1) ; for j=1:size(alleeg,2) ; for k=1:size(alleeg,3) ; smtheeg(i,j,k,:) = imfilter(squeeze(alleeg(i,j,k,:)),fspecial('gaussian',[9,1],5)) ; end ; end ; end
smtheeg = squeeze(mean(smtheeg,2)) ; 
postelecs = [23,56,24,57,25,58,26,59,27,60,61,62,63,64,29,30,31] ;

% save data: 
cd c:/shared/savedata_struct ; 
large_mersp = allamersp ; save('large_mersp','large_mersp') ; 
large_inots = inots ; save('large_inots','large_inots') ; 
large_smtheeg = smtheeg ; save('large_smtheeg','large_smtheeg') ; 
large_sumcorrs = sumcorrs ; save('large_sumcorrs','large_sumcorrs') ; 
large_dists = alldists ; save('large_dists','large_dists') ; 
large_percs = allmstims ; save('large_percs','large_percs') ; 
large_cortsz = sumgms ; save('large_cortsz','large_cortsz'); 
large_binots = binots ; save('large_binots','large_binots'); 

% compute correlation matrices: 
for i=1:30 ; for j=1:128 ; [inotcorrs(i,j),inotps(i,j)] = corr(squeeze(inots(:,i)),squeeze(mean(smtheeg(:,postelecs,j),2))) ; end ; end
for i=1:128 ; [distcorrs(i),distps(i)] = corr(mean(alldists(:,postelecs),2),mean(smtheeg(:,postelecs,i),2)) ; end 
for i=1:64 ; for j=1:128 ; [disteleccorrs(i,j),distelecps(i,j)] = corr(alldists(:,i),smtheeg(:,i,j)) ; end ; end
for i=1:64 ; [gammaeleccorrs(i),gammaelecps(i)] = corr(alldists(:,i),mean(smtheeg(:,i,40:90),3)) ; end
% compute single correlations: 
[single_inotcorr,single_inotp] = corr(mean(inots(:,8:15),2),mean(mean(smtheeg(:,postelecs,40:90),2),3)) ; 
[single_distcorr,single_distp] = corr(mean(alldists(:,postelecs),2),mean(mean(smtheeg(:,postelecs,40:90),2),3)) ; 

%peak frequency correlations 
for i=1:20 ; [peakc(i),peakp(i)] = corr(sumcorrs(peaksubs,i),allpeaks(peaksubs)') ; end
for i=1:20 ; [sumc(i),sump(i)] = corr(sumcorrs(:,i),squeeze(mean(mean(smtheeg(:,postelecs,40:90),2),3))) ; end
[peakpowr,peakpowp] = corr(allpeaks(peaksubs)',squeeze(mean(mean(smtheeg(peaksubs,postelecs,40:90),2),3))) ; 
[peakdistr,peakdistp] = corr(allpeaks(peaksubs)',mean(alldists(peaksubs,postelecs),2)) ; 
for i=1:30 ; [peakinotr(i),peakinotp(i)] = corr(allpeaks(peaksubs)',inots(peaksubs,i)) ; end

% inter-stimulus correlations:
mgamma = squeeze(mean(mean(allamersp(:,:,20:end,50:160),3),4)) ;
stimcorrs = corr(mgamma) ;

% visualization ; 
f=figure,
subplot(1,2,1) ; topoplot(gammaeleccorrs,baseeg.chanlocs,'plotchans',1:64,'numcontour',4,'conv','on','maplimits',[-.75,.75]) ; 
subplot(1,2,2) ; topoplot(squeeze(mean(mean(smtheeg(:,:,40:90),1),3)),baseeg.chanlocs,'plotchans',1:64,'numcontour',4,...
    'style','both','conv','on','emarker2',{postelecs,'o','w',5,1},'maplimits',[0,1]) ; 
figure,subplot(2,2,1) ; 
[rho,pval] = corr(mean(alldists(:,postelecs),2),mean(mean(smtheeg(:,postelecs,40:100),2),3)) ; 
plot(mean(alldists(:,postelecs),2),mean(mean(smtheeg(:,postelecs,40:100),2),3),'ok','LineWidth',2) ; lsline ; xlim([50,63]) ; ylim([-.75,2.25]) ; 
xlabel('distance(mm)') ; ylabel('gamma power (task-rest)') ; 
title(['rho=',num2str(rho),', p=',num2str(pval)]) ; 
subplot(2,2,3) ; 
[rho,pval] = corr(mean(alldists(:,postelecs),1)',squeeze(mean(mean(smtheeg(:,postelecs,40:90),1),3))') ; 
errorbarxy(mean(alldists(:,postelecs),1),squeeze(mean(mean(smtheeg(:,postelecs,40:90),1),3)),...
        std(alldists(:,postelecs),0,1)/sqrt(24),squeeze(std(mean(smtheeg(:,postelecs,40:90),3),0,1))./sqrt(24),{'.','k','k'}) ; lsline; xlim([25,80]); ylim([0.2,1.4]) ;
xlabel('distance(mm)') ; ylabel('gamma power (task-rest)') ; 
title(['rho=',num2str(rho),', p=',num2str(pval)]) ; 
subplot(2,2,4) ; 
[rho,pval] = corr(mean(alldists(:,goods),1)',squeeze(mean(mean(smtheeg(:,goods,40:90),1),3))') ; 
errorbarxy(mean(alldists(:,goods),1),squeeze(mean(mean(smtheeg(:,goods,40:90),1),3)),...
        std(alldists(:,goods),0,1)/sqrt(24),squeeze(std(mean(smtheeg(:,goods,40:90),3),0,1))./sqrt(24),{'.','k','k'}) ; lsline ;
xlabel('distance(mm)') ; ylabel('gamma power (task-rest)') ; 
title(['rho=',num2str(rho),', p=',num2str(pval)]) ; 
subplot(2,2,2) ; 
[rho,pval] = corr(mean(inots(:,15),2),mean(mean(smtheeg(:,postelecs,40:100),2),3)) ; 
plot(mean(inots(:,15),2),mean(mean(smtheeg(:,postelecs,40:100),2),3),'ok','LineWidth',2) ; lsline ; 
title(['rho=',num2str(rho),', p=',num2str(pval)]) ; ylim([-.6,2]) ; xlim([0.9,1])

figure
subplot(2,2,1) ; 
plot(sumcorrs(peaksubs,8),allpeaks(peaksubs),'ok','LineWidth',2) ; lsline ; title(['rho=',num2str(peakc(8)),', p=',num2str(peakp(8))]) ; 
xlim([-5000,35000]) ; xlabel('size of activation (#voxels)') ; ylabel('peak gamma frequency(hz)') ; 
subplot(2,2,2) ; 
plot(sumcorrs(:,8),squeeze(mean(mean(smtheeg(:,postelecs,40:90),2),3)),'ok','LineWidth',2) ; lsline ; title(['rho=',num2str(sumc(8)),', p=',num2str(sump(8))]) ; 
xlim([-5000,35000]) ; ylim([-1,2.5]); xlabel('size of activation (#voxels)') ; ylabel('gamma power (task-rest)') ; 
subplot(2,2,3) ; 
plot(squeeze(mean(mean(smtheeg(peaksubs,postelecs,40:90),2),3)),allpeaks(peaksubs)','ok','LineWidth',2) ;  lsline ; title(['rho=',num2str(peakpowr),', p=',num2str(peakpowp)]) ; 
xlabel('gamma power (task-rest)') ; ylabel('peak frequency(hz)') ; xlim([-1,2.5]) ; 
subplot(2,2,4) ;
plot(mean(alldists(peaksubs,postelecs),2),allpeaks(peaksubs)','ok','LineWidth',2) ; lsline ; title(['rho=',num2str(peakdistr),', p=',num2str(peakdistp)]) ; 
xlabel('distance(mm)') ; ylabel('peak frequency(hz)') ; xlim([48,66]) ;

figure,
subplot(2,2,1) ; shadedErrorBar(1:2:120,squeeze(mean(mean(mean(allamersp(:,:,:,60:160),1),2),4)),squeeze(std(mean(mean(allamersp(:,:,:,60:160),2),4),0,1))/sqrt(24)) ; 
hline(0,'k') ; ylabel('gamma power (db)') ; xlim([0,120]) ; xlabel('frequency(hz)') ; 
subplot(2,2,3) ; plot(distcorrs,'k','LineWidth',2) ; hline(0,'k') ; ylim([-.6,.4]) ; xlim([0,128]) ;
for i=1:length(distps) ; if distps(i) < 0.05 ; text(i,0.32,'*','Color','m') ; end ; end ; 
for i=1:length(distps) ; if distps(i) < 0.15 ; text(i,0.34,'*','Color','r') ; end ; end ; 
xlabel('frequency(hz)') ; ylabel('correlation (rho)') ; 
subplot(2,2,4) ; plot(inotcorrs(15,:),'k','LineWidth',2) ; hline(0 ,'k') ; ylim([-.6,.4]) ; xlim([0,128]) ;
for i=1:length(distps) ; if inotps(15,i) < 0.05 ; text(i,0.32,'*','Color','m') ; end ; end ; 
for i=1:length(distps) ; if inotps(15,i) < 0.15 ; text(i,0.34,'*','Color','r') ; end ; end ; 
xlabel('frequency(hz)') ; ylabel('correlation (rho)') ; 

figure,
subplot(2,2,1) ; imagesc(1:6,1:6,stimcorrs,[-1,1]) ; 
[rho13,p13] = corr(mgamma(:,1),mgamma(:,3)) ; [rho26,p26] = corr(mgamma(:,2),mgamma(:,6)) ;
subplot(2,2,3) ; plot(mgamma(:,1),mgamma(:,3),'ok','LineWidth',2) ; lsline ; title(['rho=',num2str(rho13),', p=',num2str(p13)]) ;
xlabel('stimulus 1 (unperturbed)') ; ylabel('stimulus 3 (33% contrast)') ; 
subplot(2,2,4) ; plot(mgamma(:,2),mgamma(:,6),'ok','LineWidth',2) ; lsline ; title(['rho=',num2str(rho26),', p=',num2str(p26)]) ;
xlabel('stimulus 2 (5% contrast)') ; ylabel('stimulus 6 (60% randomized)') ; 

figure,for i=1:6 ; subplot(2,3,i) ; imagesc(times,freqs,squeeze(mean(allamersp(:,i,:,:),1)),[-1.5,1.5]) ; axis xy ; if i==1 ; xlabel('time(s)') ; ylabel('frequency(hz)') ; end ; end
figure,subplot(2,2,1) ; imagesc(times,freqs,squeeze(mean(mean(allamersp(:,:,:,:),1),2)),[-1,1]) ; xlabel('time(s)') ; axis xy ; ylabel('frequency(hz)') ; 
figure,subplot(1,3,1) ; imagesc(times,freqs,squeeze(mean(allamersp(17,[1,5],:,:),2)),[-3,3]) ; axis xy ;  xlabel('time(s)') ; ylabel('frequency(hz)') ; 
subplot(1,3,2) ; imagesc(times,freqs,squeeze(mean(allamersp(4,[1,5],:,:),2)),[-3,3]) ;axis xy ; 
subplot(1,3,3) ; imagesc(times,freqs,squeeze(mean(allamersp(9,[1,5],:,:),2)),[-3,3]) ;axis xy ; 

