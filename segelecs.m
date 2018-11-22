%%% script to segment electrodes%%% process and hand label UTE images.
clear all ; close all ; 
cd c:/shared/lastute/ ;
mongs=dir('*') ; mongs(1:2) = [] ; 

for m=1:length(mongs) ; 
    
cd(['C:\shared\lastute\',mongs(m).name]) ; ls ; 
disp('loading raw data...') ; 
rute = load_untouch_nii('res_ute.nii.gz') ; ruteorig = double(rute.img) ;  
fmask = load_untouch_nii('finalmask.nii.gz') ; largemask = imdilate(fmask.img>0,strel(ones(7,7,7))) ; 
maskimg = imdilate(fmask.img>0,strel(ones(3,3,3))) ; maskimg = ~maskimg ; maskimg = maskimg.*largemask ; 
outerint = maskimg.*ruteorig ; 

boxl = 10 ; 
outerint = pad3d(outerint,boxl) ; 
boximg = zeros(size(outerint)) ; 
squareimg = zeros(size(outerint)) ; 
pointimg = zeros(size(outerint)) ; 
segimg = zeros(size(outerint)) ; 
segcoordimg = zeros(size(outerint)) ; 
rndcoordimg = zeros(size(outerint)) ; 
segcoordnumbers = zeros(size(outerint)) ; 

for c=1:30 ; 
coords = load(['mricoords_',num2str(c),'.mat']) ; 
coords = coords.mricoords + boxl ; 
for i=1:size(coords,2)
    boxi = outerint(coords(1,i)-boxl:coords(1,i)+boxl,coords(2,i)-boxl:coords(2,i)+boxl,coords(3,i)-boxl:coords(3,i)+boxl) ; 
    if c==1
        squareimg(coords(1,i)-boxl:coords(1,i)+boxl,coords(2,i)-boxl:coords(2,i)+boxl,coords(3,i)-boxl:coords(3,i)+boxl) = 1 ; 
        pointimg(coords(1,i),coords(2,i),coords(3,i)) = 1 ; 
        origcoords = coords ; 
    end
    [cx,cy,cz] = centmass3(boxi) ; normcent(c,i,:) = [cx+coords(1,i),cy+coords(2,i),cz+coords(3,i)] ;
    resboxi = reshape(boxi,[1,numel(boxi)]) ; [sv,si] = sort(resboxi,'descend') ; 
    resboxi(si(1:200)) = 1  ; resboxi(si(101:end)) = 0 ; 
    boxi = reshape(resboxi,size(boxi)) ; 
    [cx2,cy2,cz2] = centmass3(boxi) ; bincent(c,i,:) = [cx2+coords(1,i),cy2+coords(2,i),cz2+coords(3,i)] ; 
    allboxes(i,:,:,:) = boxi ; 
    if c==1 ; 
        segimg(coords(1,i)-boxl:coords(1,i)+boxl,coords(2,i)-boxl:coords(2,i)+boxl,coords(3,i)-boxl:coords(3,i)+boxl) = boxi ; 
        segcoordimg(cx2+coords(1,i)-boxl,cy2+coords(2,i)-boxl,cz2+coords(3,i)-boxl) = 1 ; 
    end
    rndcoordimg(coords(1,i),coords(2,i),coords(3,i)) = 1 ; 
end
end

stdnorms = squeeze(std(normcent,0,1)) ; 
stdbins = squeeze(std(bincent,0,1)) ; 
subnorms(m,:) = (mean(stdnorms,2)) ;
subbins(m,:) = (mean(stdbins,2)) ;

mbincent = squeeze(round(mean(bincent,1))) ;
mstdbins = mean(stdbins,2) ; 
for i=1:65 ; 
    if mstdbins(i) > 1 
        segcoordnumbers(origcoords(1,i),origcoords(2,i),origcoords(3,i))  = i ; 
    else
        segcoordnumbers(mbincent(i,1)-boxl,mbincent(i,2)-boxl,mbincent(i,3)-boxl)  = i ; 
    end
    
end

% create an image with all the centers of mass
boximg2 = zeros(size(boximg)) ; 
for i=1:size(normcent,1)
    for j=1:size(normcent,2)
        boximg(normcent(i,j,1)-boxl,normcent(i,j,2)-boxl,normcent(i,j,3)-boxl) = j ; 
        boximg2(bincent(i,j,1)-boxl,bincent(i,j,2)-boxl,bincent(i,j,3)-boxl) = j ; 
    end
end

boximg = boximg(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 
boximg2 = boximg2(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 
squareimg = squareimg(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 
pointimg = pointimg(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 
segimg = segimg(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 
segcoordimg = segcoordimg(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 
rndcoordimg = rndcoordimg(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 
segcoordnumbers = segcoordnumbers(boxl:end-boxl-1,boxl:end-boxl-1,boxl:end-boxl-1) ; 

rute.img = boximg2>0 ; save_untouch_nii(rute,'boximg2.nii.gz') ; 
rute.img = boximg>0 ; save_untouch_nii(rute,'boximg.nii.gz') ; 
rute.img = squareimg ; save_untouch_nii(rute,'squareimg.nii.gz') ; 
rute.img = imdilate(pointimg,strel(ones(3,3,3))) ; save_untouch_nii(rute,'pointimg.nii.gz') ; 
rute.img = segimg ; save_untouch_nii(rute,'segimg.nii.gz') ; 
rute.img = imdilate(segcoordimg,strel(ones(3,3,3))) ; save_untouch_nii(rute,'segcoordimg.nii.gz') ; 
rute.img = rndcoordimg ; save_untouch_nii(rute,'rndcoordimg.nii.gz') ; 
rute.img = segcoordnumbers ; save_untouch_nii(rute,'segcoordnumbers.nii.gz') ; 

end






%{

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
msubbins = mean(subbins,1) ; 
props = subbins<0.5 ; 
mprops = mean(props,1) ; 
for i=1:length(revinds)
    if revinds(i) ~= 0 
        topogms(i) = msubbins(revinds(i)) ;    
        topoprops(i) = mprops(revinds(i)) ;    
    end
end
topoplot(topogms,EEG.chanlocs,'maplimits',[0,2],'electrodes','labels')
topoplot(topoprops,EEG.chanlocs,'maplimits',[0,1],'electrodes','labels')

%}




















