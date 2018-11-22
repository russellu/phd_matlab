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

for sb=1:length(subs) ; disp(sb) ; 
    cd(['c:/shared/gamma_t1s/sub_',subs{sb},'/mri']) ;  
    locs = load('locs') ; locs = locs.locs ; 
    meancorrs = load_untouch_nii('meancorrs.nii.gz') ; 
    meancorrs.img(:,:,110:end) = 0 ; 
    [sv,si] = sort(meancorrs.img(:),'descend') ; 
    zcorrs = zeros(size(meancorrs.img)) ; zcorrs(si(1:10000)) = 1 ; 
    
    leftnorms = load_untouch_nii('surf.nii.gz') ; 
    rightnorms = load_untouch_nii('surfr.nii.gz') ; 
    bothnorms = leftnorms.img + rightnorms.img ; 
    
    [nx,ny,nz] = ind2sub(size(zcorrs),si(1:10000)) ; 
    allnorms = zeros(length(nx),3) ; 
    for i=1:length(nx)
        allnorms(i,:) = squeeze(bothnorms(nx(i),ny(i),nz(i),:)) ; 
    end
    inots(sb) = 1 - sum(abs(sum(allnorms,1)))./length(nx) ; 
    eeg = load(['c:/shared/allres/',subs{sb},'/afreps.mat']) ; eeg = eeg.afreps	 ; 
    inotmat = zeros(1,64) ; 
  %  bads = find(distmat==0) ; goods = find(distmat~=0) ; 
    alleeg(sb,:,:,:) = squeeze(mean(eeg(:,:,:,500:1000),4)) ; 
 
end

freqs = 1:4:120 ; 

clear fcorrs fps ; 
for i=1:64 ;
    for j=1:30
        [corrs(i,j),ps(i,j)] = corr(squeeze(mean(alleeg(:,:,i,j),2)),inots') ; 
    end
    [fcorrs(i),fps(i)] = corr(squeeze(mean(mean(alleeg(:,:,i,freqs>40 & freqs<100),2),4)),inots') ; 
end

subplot(1,2,1) ; topoplot(fcorrs,baseeg.chanlocs,'maplimits',[-.7,.7],'electrodes','numbers') ; 
subplot(1,2,2) ; topoplot(fps,baseeg.chanlocs,'maplimits',[0,.05],'electrodes','numbers') ; 

plot(inots,squeeze(mean(mean(alleeg(:,:,23,freqs>40 & freqs<100),2),4)),'o') ; lsline ; 
[rho,p] = corr(inots',squeeze(mean(mean(alleeg(:,:,23,freqs>40 & freqs<100),2),4))) ; 
title(['rho = ',num2str(rho),', p=',num2str(p)]) ; xlabel('i-not') ; ylabel('gamma power') ; 











