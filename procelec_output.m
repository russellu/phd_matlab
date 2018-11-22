clear all ; close all ; 
cd c:/shared/lastute/ ;
mongs=dir('*') ; mongs(1:2) = [] ; 

csize = 50 ; 
for mong=1:length(mongs) ; 
    cd(['C:\shared\lastute\',mongs(mong).name]) ; ls ; 
    padamt = 15 ; 
    for p=1:20 ; disp(p) ; 
        orignii = load_untouch_nii(['origimg',num2str(p),'.nii.gz']) ; 
        coords = load(['mricoords_',num2str(p)]) ; coords = coords.mricoords ; 
        allcoords(p,:,:) = coords ; 
        img = orignii.img ; img = pad3d(img,padamt) ; 
        coords = coords + padamt  ;
        esize = 8  ;
        newimg = zeros(size(img)) ; 
        for i=1:size(coords,2)
            boxes(i,:,:,:) = img(coords(1,i)-esize:coords(1,i)+esize,coords(2,i)-esize:coords(2,i)+esize,coords(3,i)-esize:coords(3,i)+esize) ;  
            boxi = squeeze(boxes(i,:,:,:)) ; 
            [sv,si] = sort(boxi(:),'descend') ; 
            boxi = zeros(size(boxi)) ; 
            boxi(si(1:csize)) = 1 ; 
            newimg(coords(1,i)-esize:coords(1,i)+esize,coords(2,i)-esize:coords(2,i)+esize,coords(3,i)-esize:coords(3,i)+esize) = boxi*i ; 
        end
        orig = newimg(padamt:end-(padamt+1),padamt:end-(padamt+1),padamt:end-(padamt+1)) ; 
        if p==1 ; orig1 = orig ; end
        for i=1:length(coords)
           [cx,cy,cz] = centmass3(orig==i) ;  
            newcoords(i,p,:) = [cx,cy,cz] ; 
        end
    end
    figure,
    plot(sort(squeeze(mean(std(newcoords,0,2),3)))) ; title(num2str(mean(mean(std(newcoords,0,2))))) ; hold on ; 
    plot(sort(squeeze(mean(std(allcoords,0,1),2))),'r') ; 
    
    subnew(mong,:,:,:) = newcoords ; 
    suball(mong,:,:,:) = allcoords ; 
    
    %subnewsizes(csizecount,mong,:,:,:) = newcoords ; 
    coordbrain = zeros(size(orignii.img)) ; 
    segcoordbrain = zeros(size(orignii.img)) ; 
    labcoordbrain = zeros(size(orignii.img)) ; 
    
    coords = squeeze(allcoords(1,:,:)) ; 
    segcoords = squeeze(newcoords(:,1,:))' ; 
    for i=1:size(coords,2) ; 
        coordbrain(coords(1,i),coords(2,i),coords(3,i)) = 1 ;
        segcoordbrain(segcoords(1,i),segcoords(2,i),segcoords(3,i)) = 1 ;
        labcoordbrain(segcoords(1,i),segcoords(2,i),segcoords(3,i)) = i ;
    end
    coordbrain  = imdilate(coordbrain,strel(ones(3,3,3))) ; 
    segcoordbrain  = imdilate(segcoordbrain,strel(ones(3,3,3))) ; 

    orignii.img = coordbrain ; save_untouch_nii(orignii,'handcoordbrain.nii.gz') ; 
    orignii.img = segcoordbrain ; save_untouch_nii(orignii,'segcoordbrain.nii.gz') ; 

    orignii.img = orig1>0 ; save_untouch_nii(orignii,'segs.nii.gz') ; 

end

%{
%csizecount = csizecount + 1 ; 
%end
stdsubs = squeeze(std(subnew,0,3)) ; 
mstd = squeeze(mean(mean(stdsubs,1),3)) ; 
mdirs = squeeze(mean(stdsubs,3)) ; 

EEG = pop_loadset('C:\shared\badger_eeg\alex\1hz_preproc_retino_allstims_01_Pulse Artifact Correction.set') ;
elabs = {EEG.chanlocs.labels} ;
elecorder = {'FP1','FPZ','FP2','AF8','AF4','GND','AF3','AF7','F7','F5','F3','F1','FZ','F2','F4','F6','F8','FT10','FT8','FC6','FC4','FC2','REF','FC1','FC3','FC5','FT7','FT9',...
    'T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP10','TP8','CP6','CP4','CP2','CPZ','CP1','CP3','CP5','TP7','TP9','P7','P5','P3','P1','PZ','P2','P4','P6','P8',...
    'PO8','PO4','POZ','PO3','PO7','O1','OZ','O2'} ; 

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


propdirs = mdirs<.5 ; 
props = sum(propdirs,1)./14 ; 


mrand = squeeze(mean(mean(std(suball,0,2),1),3)) ; 


clear topogms topopercs topornd
for i=1:length(revinds)
    if revinds(i) ~= 0 
        topogms(i) = mstd(revinds(i)) ;    
        topopercs(i) = props(revinds(i)) ;
        topornd(i) = mrand(revinds(i)) ; 
    end
end

topoplot(topogms,EEG.chanlocs,'maplimits',[0,2],'electrodes','labels')

topoplot(topopercs,EEG.chanlocs,'maplimits',[0,1],'electrodes','labels')

topoplot(topornd,EEG.chanlocs,'maplimits',[0,2],'electrodes','labels')



barh([topogms;topornd]') ; vline(1,'k') ; set(gca,'YTick',1:64,'YTickLabel',elabs) ; 
%}
