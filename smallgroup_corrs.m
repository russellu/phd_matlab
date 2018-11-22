clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'} ; 
elecorder = load('E:\all_white_normals\a2_good\elecorder') ; elecorder = elecorder.elecorder ; 
baseeg = pop_loadset('E:\badger_eeg\alex\retino_rest.set') ;
peaksubs = [1,2,3,4,5,6,8,9] ; 
allpeaks = [50,54,58,60,56,52,64,54,56] ; %64
elabs = {baseeg.chanlocs.labels} ; 
for sb=1:length(subs) ; disp(sb) ; 
    cd(['E:/badger_struct/',subs{sb}]) ; 
    %mpercs = load('mpercs') ; mpercs = mpercs.mpercs; allmpercs(sb,:,:) = mpercs; 
    gm = load_untouch_nii('gm.nii.gz'); 
    sumgms(sb) = sum(gm.img(:)); 
    meants = load(['E:/badger_struct/',subs{sb},'/meants']) ; meants = meants.meants ; allmeants(sb,:,:) = meants ;  
    mstersp = load('mstersp') ; mstersp = mstersp.mstersp ; allamersp(sb,:,:,:) = mstersp ; 
    corrbrain = load_untouch_nii('meancorrs_fs.nii.gz') ; 
    threshs = 0.1:0.05:0.9 ; 
    for thr=1:length(threshs)
        threshcorrs = corrbrain.img>threshs(thr) ; sumcorrs(sb,thr) = sum(threshcorrs(:)) ; 
    end
    t1 = load_untouch_nii('T1.nii.gz') ; 
    [sv,si] = sort(corrbrain.img(:),'descend') ; 
    zimg = zeros(size(corrbrain.img)) ; zimg(si(1:8000)) = 1 ; 
    zimg(:,1:125,:) = 0 ; 
    [cx,cy,cz] = centmass3(zimg) ; 
    [bsv,bsi] = sort(zimg(:),'descend') ; 
    [sx,sy,sz] = ind2sub(size(corrbrain.img),bsi(1:8000)) ; 
    coords = load('coords') ; coords = coords.coords ; 
    sqrdiffs = sqrt((cx-coords(:,1)).^2 + (cy-coords(:,2)).^2 + (cz-coords(:,3)).^2) ; 
    for i=1:size(coords,1)
        for j=1:length(sx)
            mindiffs(i,j) = sqrt((sx(j)-coords(i,1)).^2 + (sy(j)-coords(i,2)).^2 + (sz(j)-coords(i,3)).^2) ; 
        end
    end
    mindiffs = min(mindiffs,[],2) ;  
    dists = zeros(size(elabs)) ; 
    for i=1:length(sqrdiffs) ; if ~isempty(find(strcmpi(elecorder{i},elabs))) ; elecind = find(strcmpi(elecorder{i},elabs)) ; dists(elecind) = sqrdiffs(i) ; end ; end
    alldists(sb,:) = dists ; 
    
    %inot
    leftnorms = load_untouch_nii('left_norms.nii.gz') ; 
    rightnorms = load_untouch_nii('right_norms.nii.gz') ; 
    bothnorms = leftnorms.img + rightnorms.img ; 
    normls = 100:500:15000 ; 
    for nl=1:length(normls) 
    [nx,ny,nz] = ind2sub(size(zimg),si(1:normls(nl))) ; 
    allnorms = zeros(length(nx),3) ; 
    for i=1:length(nx)
        allnorms(i,:) = squeeze(bothnorms(nx(i),ny(i),nz(i),:)) ; 
    end
    inots(sb,nl) = 1 - sum(abs(sum(allnorms,1)))./length(nx) ; 
    end  
end


cd E:\badger_eeg\valerie\outside
mampdiff = load('mampdiff'); mampdiff = mampdiff.mampdiff; 
postchans = [27,51,7,45];
for i=1:64
   [ecorrs(i),eps(i)] = corr(squeeze(alldists(:,i)),squeeze(mean(mampdiff(:,i,20:end),3)));  
    
end

subplot(2,2,1)
plot(squeeze(mean(alldists(:,postchans),2)),squeeze(mean(mean(mampdiff(:,postchans,20:end),2),3)),'kd','LineWidth',2); lsline
[c,p] = corr(squeeze(mean(alldists(:,postchans),2)),squeeze(mean(mean(mampdiff(:,postchans,20:end),2),3))); 
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]); xlabel('distance(mm)'); ylabel('gamma amplitude'); 
subplot(2,2,2); 
topoplot(ecorrs,baseeg.chanlocs,'maplimits',[-.7,.7],'style','both','emarker2',{find(eps<0.05),'.','w'}); colormap parula;

subplot(2,2,3); 
%errorbarxy(mean(alldists,1),squeeze(mean(mean(mampdiff(:,:,20:end),1),3)),std(alldists,0,1)./sqrt(9),squeeze(std(mean(mampdiff(:,:,20:end),3),0,1))/3,{'kd'}); 






