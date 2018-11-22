clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 
elabs = load('C:\shared\all_white_normals\a2_good\elabs') ; elabs = elabs.elabs ; 
elecorder = load('C:\shared\all_white_normals\a2_good\elecorder') ; elecorder = elecorder.elecorder ; 
baseeg = pop_loadset('c:/shared/allres/alex/merged.set') ; 
baselabs = {baseeg.chanlocs.labels} ;
for i=1:length(baselabs) ;
    if ~isempty(find(strcmpi(baselabs{i},elecorder)))
        orderinds(i) = find(strcmpi(baselabs{i},elecorder)) ;
    end
end 
% first - make a topoplot of distance from V1 for each subject 
% orderinds - the location index of each eeg structure electrode

for sb=1:length(subs)
    cd(['c:/shared/gamma_t1s/sub_',subs{sb},'/mri']) ;  
    locs = load('locs') ; locs = locs.locs ; 
    meancorrs = load_untouch_nii('meancorrs.nii.gz') ;     t1 = load_untouch_nii('brain.nii.gz') ; 

    allmeancorrs(:,:,:,sb) = (meancorrs.img>.2).*(t1.img>0) ; 
    meancorrs.img(:,:,110:end) = 0 ; 
    allt1s(:,:,:,sb) = t1.img ; 
    
    nvoxes(sb) = sum(meancorrs.img(:)>0.5) ;
    
    %subplot(3,8,sb) ; imagesc(squeeze(max(meancorrs.img,[],2))) ; 
    [sv,si] = sort(meancorrs.img(:),'descend') ; 
    zcorrs = zeros(size(meancorrs.img)) ; zcorrs(si(1:10000)) = 1 ; 
    [cx,cy,cz] = centmass3(zcorrs) ; 
    sqrdiffs = sqrt((locs(1,:)-cx).^2 + (locs(2,:)-cy).^2 + (locs(3,:)-cz).^2) ; 
    
    eeg = load(['c:/shared/allres/',subs{sb},'/afreps.mat']) ; eeg = eeg.afreps	 ; 
    distmat = zeros(1,64) ; 
    for i=1:length(orderinds) ; if orderinds(i)~=0 ; distmat(i) = sqrdiffs(orderinds(i)) ; end ; end
    bads = find(distmat==0) ; goods = find(distmat~=0) ; 
    alldistmat(sb,:) = distmat ; 
    alleeg(sb,:,:,:) = squeeze(mean(eeg(:,:,:,500:1000),4)) ; 
    allteeg(sb,:,:,:,:) = eeg ;  
    allsqrdiffs(sb,:) = sqrdiffs ; 
    %subplot(4,8,sb) ; topoplot(distmat,baseeg.chanlocs,'maplimits',[0,150],'plotchans',goods) ; 
end

for i=1:size(alleeg,2)
    for j=1:size(alleeg,3)
        for k=1:size(alleeg,4)
            corrs(i,j,k) = corr2(squeeze(alleeg(:,i,j,k)),squeeze(alldistmat(:,j))) ; 
        end
    end
end
clear corrs ; 
for i=1:64
    for j=1:30
        corrs(i,j) = corr2(squeeze(mean(alleeg(:,:,i,j),2)),squeeze(alldistmat(:,i))) ; 
    end
end
%for i=1:30 ; subplot(5,6,i) ; topoplot(squeeze((corrs(:,i))),baseeg.chanlocs,'maplimits',[-.5,.5],'plotchans',goods) ; end
freqs = 1:4:120 ; 
for i=1:64 ; %subplot(5,13,i) ; 
 %   plot(squeeze(mean(mean(alleeg(:,:,i,15:25),2),4)),alldistmat(:,i),'.') ; lsline ;
    [bcorrs(i),bps(i)] = corr(squeeze(mean(mean(alleeg(:,:,i,freqs>40 & freqs<100),2),4)),alldistmat(:,i)) ; %title(num2str(bcorrs(i))) ; 
end
subplot(1,2,1) ; topoplot(squeeze((bcorrs)),baseeg.chanlocs,'maplimits',[-.7,.7],'plotchans',goods,'electrodes','numbers') ; es = [24,57,61] ; 
subplot(1,2,2) ; topoplot(squeeze((bps)),baseeg.chanlocs,'maplimits',[0,.05],'plotchans',goods,'electrodes','labels') ;

plot(mean(alldistmat(:,es),2),squeeze(mean(mean(mean(alleeg(:,:,es,freqs>40 & freqs<100),2),4),3)),'o') ; lsline ; 
[rho,p] = corr(squeeze(mean(mean(mean(alleeg(:,:,es,freqs>40 & freqs<100),2),4),3)),mean(alldistmat(:,es),2)); 
title(['rho = ',num2str(rho),', p=',num2str(p)]) ; xlabel('distance(mm)') ; ylabel('gamma power') ; 

for i=1:size(alleeg,3)
    for j=1:size(alleeg,4)
        [vrho(i,j),vp(i,j)]= corr(squeeze(mean(alleeg(:,:,i,j),2)),nvoxes') ; 
    end
end

subplot(1,2,1) ; topoplot(squeeze(mean(vrho(:,freqs>40 & freqs<100),2)),baseeg.chanlocs,'maplimits',[-.7,.7],'plotchans',goods,'electrodes','labels') ;
subplot(1,2,2) ; topoplot(squeeze(mean(vrho(:,freqs>40 & freqs<100),2)),baseeg.chanlocs,'maplimits',[0,0.05],'plotchans',goods,'electrodes','labels') ;

plot(nvoxes,squeeze(mean(mean(mean(alleeg(:,:,[61,62,63,29,39,31],freqs>40 & freqs<100),2),4),3)),'o') ; xlim([-1000,35000]) ; lsline ; 
[rho,p] = corr(squeeze(mean(mean(mean(alleeg(:,:,[61,62,63,29,39,31],freqs>40 & freqs<100),2),4),3)),nvoxes') ; 
title(['rho = ',num2str(rho),', ','p=',num2str(p)]) ; 
ylabel('gamma power') ; xlabel('number of voxels') ; 

figure,
for i=1:24;  subplottight(2,12,i) ; 
    allmeancorrs(1,80:200,1,i) = 0 ; 
    plotoverlayIntensity2D(squeeze(allt1s(:,150,:,i)),squeeze((mat2gray(mean(allmeancorrs(:,80:200,:,i),2)))),squeeze(mean(allmeancorrs(:,80:200,:,i),2)),90) ; 
end

