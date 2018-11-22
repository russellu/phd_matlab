clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'} ; 
%elabs = load('c:/shared/elabs') ; elabs = elabs.elabs ; 
elecorder = load('c:/shared/elecorder') ; elecorder = elecorder.elecorder ; 
baseeg = pop_loadset(['c:/shared/simdenoise/alex','/bsclean_newgrad_retino_gamma_0',num2str(1),'.set']) ;
elabs = {baseeg.chanlocs.labels} ; 
for sb=1:length(subs) ; disp(sb) ; 
    cd(['c:/shared/newbadger_mri/',subs{sb}]) ; 
    meants = load(['c:/shared/simdenoise/',subs{sb},'/meants']) ; meants = meants.meants ; allmeants(sb,:,:) = meants ;  
    corrbrain = load_untouch_nii('trimcorrs.nii.gz') ; 
    [sv,si] = sort(corrbrain.img(:),'descend') ; 
    zimg = zeros(size(corrbrain.img)) ; zimg(si(1:8000)) = 1 ; % zimg = medfilt3(zimg) ;
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
    mindists = zeros(size(elabs)) ; 
    for i=1:length(mindiffs) ; if ~isempty(find(strcmpi(elecorder{i},elabs))) ; elecind = find(strcmpi(elecorder{i},elabs)) ; mindists(elecind) = mindiffs(i) ; end ; end

    
    figure,topoplot(dists,baseeg.chanlocs,'maplimits',[0,100]) ; 
    alldists(sb,:) = dists ; allmindists(sb,:) = mindists ; 
end
smtheeg = zeros(size(allmeants)) ; 
for i=1:9 ; for j=1:64 ; smtheeg(i,j,:) = imfilter(squeeze(allmeants(i,j,:)),fspecial('gaussian',[9,1],3)) ; end ; end
for i=1:64
    for j=1:128
        corrs(i,j) = corr(alldists(:,i),smtheeg(:,i,j)) ; 
        mincorrs(i,j) = corr(allmindists(:,i),smtheeg(:,i,j)) ; 
    end
end
for i=1:64 
   [cm(i),pm(i)] = corr(alldists(:,i),mean(smtheeg(:,i,40:90),3)) ; 
   [mincm(i),minpm(i)] = corr(allmindists(:,i),mean(smtheeg(:,i,40:90),3)) ;
end

bades = [61,62,29,30] ; goodes = zeros(1,64) ; goodes(bades) = 1 ; goods = find(goodes==0) ; 
subplot(2,2,1) ; topoplot(cm(1,1:64),baseeg.chanlocs,'maplimits',[-.75,.75],'plotchans',goods) ; 
subplot(2,2,2) ; topoplot(pm(1,1:64),baseeg.chanlocs,'maplimits',[0,.05],'plotchans',goods) ; 
postelecs = [15,51,7,37,19,38,8,52,16,59,45,31,46,60,9,20,10] ; 
[c,p] = corr(mean(alldists(:,postelecs(1:end)),2),mean(mean(smtheeg(:,postelecs(1:end),40:90),2),3)) ; 
subplot(2,2,3) ; plot(mean(alldists(:,postelecs(1:end)),2),mean(mean(smtheeg(:,postelecs(1:end),40:90),2),3),'o') ; lsline ; title(['rho= ',num2str(c),', p=',num2str(p),', n=9']) ; 
ylim([-.75,2.75]) ; 
subplot(2,2,4) ; 
errorbarxy(squeeze(mean(alldists,1)),squeeze(mean(mean(smtheeg(:,:,40:90),1),3)),squeeze(std(alldists,0,1))/3,squeeze(std(mean(smtheeg(:,:,40:90),3),0,1))/3,{'o','k','b'})
xlabel('distance (mm)') ; ylabel('gamma power (task-rest)') ; 
















