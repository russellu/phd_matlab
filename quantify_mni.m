clear all ; close all ; 
EEG = pop_loadset('C:\shared\badger_eeg\alex\1hz_preproc_retino_allstims_01_Pulse Artifact Correction.set') ; 

subs = {'alex','biz','cloud','dave','felix','jeremie','karl','nic','pierre','russ','sukh','terry'} ; 
for s=1:length(subs)
    cd(['c:/shared/utef/',subs{s},'/ants_ute']) ; 
    subsi = load_untouch_nii('color_locs_in_mni.nii.gz') ; 
    imgs(s,:,:,:) = subsi.img ;
end
cd c:/shared/ATLASES ; ls 
mni = load_untouch_nii('mni_elecs.nii.gz') ; mnielecs = mni.img ; 

for i=1:size(imgs,1)
    for j=1:65 ;  
        [cx,cy,cz] = centmass3(squeeze(imgs(i,:,:,:))==j) ; 
        subelecs(i,j,:) = [cx,cy,cz] ; 
    end
end

interpsubs = zeros(size(subelecs)) ; 
for i=1:size(subelecs,2)
    for j=1:3
        elecij = squeeze(subelecs(:,i,j)) ; 
        nans = find(isnan(elecij)) ;
        if ~isempty(nans)    
            elecij(nans) = mean(elecij((~isnan(elecij)))) ; 
        end
        interpsubs(:,i,j) = elecij ;
    end
end

msubelecs = round(squeeze(mean(interpsubs))) ; 
stdsubelecs = mean(squeeze(std(interpsubs,0,1)),2) ; 

%cd c:/brainstorm3 ; labs = load('labs.mat') ; labs = labs.labs ; 
cd c:/users/butr2901/documents ; ls 
fid = fopen('coords.sxyz')  ;
data = textscan(fid,'%f %f %f %s') ; 
fclose(fid) ;
xyz = [data{1},data{2},data{3}] ; labs = data{4} ; 

cd c:/shared/utef/jeremie ; ls ; elecorder = load('elecorder.mat') ; elecorder = elecorder.elecorder ; 

for i=1:length(labs)
    [cx,cy,cz] = centmass3(mnielecs==i) ; 
    mnicoords(i,:) = [cx,cy,cz] ; 
end

clear subinds
for i=1:length(labs) ;
    if ~strcmpi(labs{i},{'FCz','F9','F10'}) ;
    subinds(i) = find(strcmpi(labs{i},elecorder)) ;
    end
end
% get the single subject x,y,z coordinates for each electrodes
for i=1:size(imgs,1)
    for j=1:length(elecorder)
        [cx,cy,cz] = centmass3(squeeze(imgs(i,:,:,:))==j) ; 
        coords(i,j,:) = [cx,cy,cz] ; 
    end
end

% subinds holds the indices of elecorder (native space) to which each MNI
% space electrode belongs
mcoords = msubelecs ; 
meanmni = zeros(size(mni.img)) ; colormni = zeros(size(mni.img)) ; 
for i=1:size(msubelecs,1)
    meanmni(msubelecs(i,1),msubelecs(i,2),msubelecs(i,3)) = 1 ;
    colormni(msubelecs(i,1),msubelecs(i,2),msubelecs(i,3)) = i ;
end
meanmni = imdilate(meanmni,strel(ones(3,3,3))) ; 
colormni = imdilate(colormni,strel(ones(3,3,3))) ; 

for i=1:length(subinds)
    if subinds(i) ~= 0
        mlocs(i,:) = mcoords(subinds(i),:) ; 
        ediffs(i) = sqrt(sum((mnicoords(i,:)-mcoords(subinds(i),:)).^2)) ;     
        stds(i) = stdsubelecs(subinds(i)) ; 
    end
end

olabs = {EEG.chanlocs.labels} ; 
for i=1:length(olabs)
    if ~strcmpi('ECG',olabs{i}) ;
        olabinds(i) = find(strcmpi(olabs{i},labs)) ; 
    end
end

stdcoords = squeeze(std(coords,0,1)) ; 

for i=1:length(olabinds)
    if olabinds(i) ~= 0 
        eeglabs(i) = ediffs(olabinds(i)) ; 
        eegstds(i) = stds(olabinds(i)) ; 
        eegstdcoords(i,:) = stdcoords(olabinds(i),:) ; 
    end
end

% get voxel vectors between each mean electrode location and the mni locations
linevol = zeros(size(mni.img)) ; 
for i=1:size(mlocs,1)
    xdiffi = mlocs(i,1) - mnicoords(i,1) ;  
    ydiffi = mlocs(i,2) - mnicoords(i,2) ;  
    zdiffi = mlocs(i,3) - mnicoords(i,3) ;  
    disti = sqrt(xdiffi.^2 + ydiffi.^2 + zdiffi.^2) ; 
    distvec = 1:ceil(disti) ; 
    uniti = [xdiffi,ydiffi,zdiffi]./disti ; 
    projvec = repmat(distvec,[3,1]).*repmat(uniti',[1,size(distvec,2)]) ;
    posi = [mlocs(i,1),mlocs(i,2),mlocs(i,3)] ; 
    projvals = repmat(posi',[1,size(projvec,2)]) - projvec ; 
    lineinds{i} = sub2ind(size(mni.img),round(projvals(1,:)),round(projvals(2,:)),round(projvals(3,:))) ; 
    linevol(lineinds{i}) = 1 ; 
end




cd c:/shared/utef/jeremie ; ls ; elecorder = load('elecorder.mat') ; elecorder = elecorder.elecorder ; 
elabs = {EEG.chanlocs.labels} ;
% you want the electrode indices of the elecorder
for i=1:length(elecorder)
   if ~strcmp(elecorder{i},'REF') && ~strcmp(elecorder{i},'GND') ; 
       elecorderinds(i) = find(strcmpi(elecorder{i},elabs)) ; 
   end
end

for i=1:length(elabs)
    if ~strcmpi(elabs{i},'ECG')
        revinds(i) = find(strcmpi(elabs{i},elecorder)) ; 
    end
end


%mnstds = squeeze(std(ndiffs,0,1)) ; mndists = squeeze(mean(ndiffs,1)) ; 
%mnstds_comp = squeeze(std(ndiffs_comp,0,1)) ; 
%msstds = squeeze(std(sdiffs,0,1)) ; msdists = squeeze(mean(sdiffs,1)) ; 
%msstds_comp = squeeze(std(sdiffs_comp,0,1)) ; 
clear topodiffs
for i=1:length(revinds)
    if revinds(i) ~= 0 
        %{
        toponstds(i,:) = mnstds(revinds(i),:) ;     
        toponstds_comp(i,:,:) = mnstds_comp(revinds(i),:,:) ; 
        
        toposstds(i,:) = msstds(revinds(i),:) ;  
        toposstds_comp(i,:,:) = msstds_comp(revinds(i),:,:) ; 
        
        topondiffs(i,:) = mndists(revinds(i),:) ; 
        toposdiffs(i,:) = msdists(revinds(i),:) ; 
        %}
        topostdcoords(i,:) = stdcoords(revinds(i),:) ; 
    end
end



subplot(3,3,7) ; topoplot(squeeze(topostdcoords(:,1)),EEG.chanlocs,'maplimits',[0,12]) ; colorbar ; title('MNI locs std(x)') ; 
subplot(3,3,8) ; topoplot(squeeze(topostdcoords(:,2)),EEG.chanlocs,'maplimits',[0,12]) ; colorbar ; title('MNI locs std(y)') ; 
subplot(3,3,9) ; topoplot(squeeze(topostdcoords(:,3)),EEG.chanlocs,'maplimits',[0,12]) ; colorbar ; title('MNI locs std(z)') ; 

subplot(2,2,1) ; barwitherr(squeeze(std(topostdcoords,0,1)),squeeze(mean(topostdcoords,1))) ; set(gca,'XTick',[1,2,3],'XTickLabel',{'x','y','z'}) ; xlabel('dimension') ; ylabel('std(mm)') ; ylim([0,10]) ; 
title('x,y,z std (MNI locations)') ; 



