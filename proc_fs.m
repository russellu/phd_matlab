%%% process and hand label UTE images.
%clear all ; close all ;
slugs = {'alex','biz','cloud','dave','felix','jeremie','karl','nic','pierre','russ','sukh','terry'} ; 
EEG = pop_loadset('C:\shared\badger_eeg\alex\1hz_preproc_retino_allstims_01_Pulse Artifact Correction.set') ; 
cd c:/shared/utef/jeremie ; ls ; elecorder = load('elecorder.mat') ; elecorder = elecorder.elecorder ; 
elabs = {EEG.chanlocs.labels} ;

for slug=1:length(slugs)
    cd(['c:/shared/utef/',slugs{slug},'/labs']) ;  
    disp(slugs{slug}) ; 
    labs=dir('*gz') ; 
    for lab=1:length(labs)
        labi = load_untouch_nii(labs(lab).name) ; 
        [cx,cy,cz] = centmass3(labi.img) ; 
        allcents(slug,lab,:) = [cx,cy,cz] ; 
    end
    cd .. ; 
    cwarps = load_untouch_nii('colorwarpelecs.nii.gz') ; % avg locations in native space
    standardcoords{slug} = cwarps.img ; 
    cnatives = load_untouch_nii('native_color_locs.nii.gz') ; % native locations in native space
    nativecoords{slug} = cnatives.img ; 
end

for i=1:length(standardcoords)
    for j=1:65 
        [cx,cy,cz] = centmass3(standardcoords{i}==j) ;
        standardlocs(i,j,:) = [cx,cy,cz] ; 
        [cx,cy,cz] = centmass3(nativecoords{i}==j) ;
        nativelocs(i,j,:) = [cx,cy,cz] ; 
    end
end
figure,imagesc(uint8(mat2gray(standardlocs)*255))
for i=1:12
    for j=1:3
        colij = squeeze(standardlocs(i,:,j)) ; 
        meansubs(i,j) = mean(colij(~isnan(colij))) ; 
    end
end
slocs2 = standardlocs ; 
for i=1:12 ;
    for j=1:64 ;
        for k=1:3 ;
            if isnan(standardlocs(i,j,k))
                rowjk = squeeze(standardlocs(:,j,k)) ; 
                meanjk = mean(rowjk(~isnan(rowjk))) ; 
                slocs2(i,j,k) = meanjk + (meansubs(i,k)-mean(meansubs(:,k))) ; 
            end
        end
    end
end
figure,imagesc(uint8(mat2gray(slocs2)*255))


for i=1:12
    for j=1:65
        for k=1:68
            sdiffs(i,j,k) = squeeze(sqrt((slocs2(i,j,1)-allcents(i,k,1)).^2 + (slocs2(i,j,2)-allcents(i,k,2)).^2 + (slocs2(i,j,3)-allcents(i,k,3)).^2)) ; % standard space difference to all FS
            sdiffs_comp(i,j,k,1) = squeeze(sqrt((slocs2(i,j,1)-allcents(i,k,1)).^2)) ;
            sdiffs_comp(i,j,k,2) = squeeze(sqrt((slocs2(i,j,2)-allcents(i,k,2)).^2)) ;
            sdiffs_comp(i,j,k,3) = squeeze(sqrt((slocs2(i,j,3)-allcents(i,k,3)).^2)) ;
           
            ndiffs(i,j,k) = squeeze(sqrt((nativelocs(i,j,1)-allcents(i,k,1)).^2 + (nativelocs(i,j,2)-allcents(i,k,2)).^2 + (nativelocs(i,j,3)-allcents(i,k,3)).^2)) ; % native space difference to all FS
            ndiffs_comp(i,j,k,1) = squeeze(sqrt((nativelocs(i,j,1)-allcents(i,k,1)).^2)) ; 
            ndiffs_comp(i,j,k,2) = squeeze(sqrt((nativelocs(i,j,2)-allcents(i,k,2)).^2)) ; 
            ndiffs_comp(i,j,k,3) = squeeze(sqrt((nativelocs(i,j,3)-allcents(i,k,3)).^2)) ; 
        end
    end
end

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


mnstds = squeeze(std(ndiffs,0,1)) ; mndists = squeeze(mean(ndiffs,1)) ; 
mnstds_comp = squeeze(std(ndiffs_comp,0,1)) ; 
msstds = squeeze(std(sdiffs,0,1)) ; msdists = squeeze(mean(sdiffs,1)) ; 
msstds_comp = squeeze(std(sdiffs_comp,0,1)) ; 
clear topodiffs
for i=1:length(revinds)
    if revinds(i) ~= 0 
        toponstds(i,:) = mnstds(revinds(i),:) ;     
        toponstds_comp(i,:,:) = mnstds_comp(revinds(i),:,:) ; 
        
        toposstds(i,:) = msstds(revinds(i),:) ;  
        toposstds_comp(i,:,:) = msstds_comp(revinds(i),:,:) ; 
        
        topondiffs(i,:) = mndists(revinds(i),:) ; 
        toposdiffs(i,:) = msdists(revinds(i),:) ; 
    end
end
clear meanstds ; 
for i=1:size(toposstds,1)
    [sv,si] = sort(toposstds(i,:),'ascend') ; 
    meannstds(i) = squeeze(mean(toponstds(i,si(1:5)))) ; 
    meannstds_comp(i,:) = squeeze(mean(toponstds_comp(i,si(1:5),:),2)) ; 
    meansstds(i) = squeeze(mean(toposstds(i,si(1:5)))) ; 
    meansstds_comp(i,:) = squeeze(mean(toposstds_comp(i,si(1:5),:),2)) ; 
end
meansstds(29) = mean(meansstds) ; meansstds_comp(29,:) = mean(meansstds_comp,1) ; 
meannstds(29) = mean(meannstds) ; meannstds_comp(29,:) = mean(meannstds_comp,1) ; 
subplot(1,3,1) ; topoplot(double(meansstds),EEG.chanlocs,'maplimits',[0,8],'electrodes','labels') ;
subplot(1,3,2) ; topoplot(double(meannstds),EEG.chanlocs,'maplimits',[0,8],'electrodes','labels') ;
%subplot(1,2,2) ; topoplot(double((meannstds-meansstds) > 0),EEG.chanlocs,'maplimits',[-1,1],'electrodes','labels') ;
subplot(1,3,3) ; topoplot(double(meannstds)-double(meansstds),EEG.chanlocs,'maplimits',[-4,4],'electrodes','labels') ;

subplot(1,3,1) ; topoplot(double(squeeze(meannstds_comp(:,1))),EEG.chanlocs,'maplimits',[0,12]) ; colorbar ; title('native locs std(x)') ; 
subplot(1,3,2) ; topoplot(double(squeeze(meannstds_comp(:,2))),EEG.chanlocs,'maplimits',[0,12]) ; colorbar ; title('native locs std(y)') ; 
subplot(1,3,3) ; topoplot(double(squeeze(meannstds_comp(:,3))),EEG.chanlocs,'maplimits',[0,12]) ; colorbar ; title('native locs std(z)') ; 
subplot(1,3,1) ; topoplot(double(squeeze(meansstds_comp(:,1))),EEG.chanlocs,'maplimits',[0,12]) ; colorbar ; title('mean locs std(x)') ;
subplot(1,3,2) ; topoplot(double(squeeze(meansstds_comp(:,2))),EEG.chanlocs,'maplimits',[0,12]) ; colorbar ; title('mean locs std(y)') ;
subplot(1,3,3) ; topoplot(double(squeeze(meansstds_comp(:,3))),EEG.chanlocs,'maplimits',[0,12]) ; colorbar ; title('mean locs std(z)') ;
subplot(1,3,1) ; topoplot(squeeze(topostdcoords(:,1)),EEG.chanlocs,'maplimits',[0,12]) ; colorbar ; title('MNI locs std(x)') ; 
subplot(1,3,2) ; topoplot(squeeze(topostdcoords(:,2)),EEG.chanlocs,'maplimits',[0,12]) ; colorbar ; title('MNI locs std(y)') ; 
subplot(1,3,3) ; topoplot(squeeze(topostdcoords(:,3)),EEG.chanlocs,'maplimits',[0,12]) ; colorbar ; title('MNI locs std(z)') ; 


figure,
subplot(2,2,1) ; barwitherr(squeeze(std(meannstds_comp,0,1)),squeeze(mean(meannstds_comp,1))) ; set(gca,'XTick',[1,2,3],'XTickLabel',{'x','y','z'}) ; xlabel('dimension') ; ylabel('std(mm)') ; ylim([0,10]) ; 
title('x,y,z std (native locations)') ; 
subplot(2,2,2) ; barwitherr(squeeze(std(meansstds_comp,0,1)),squeeze(mean(meansstds_comp,1))) ; set(gca,'XTick',[1,2,3],'XTickLabel',{'x','y','z'}) ; xlabel('dimension') ; ylabel('std(mm)') ; ylim([0,10]) ; 
title('x,y,z std (mean locations)') ; 






