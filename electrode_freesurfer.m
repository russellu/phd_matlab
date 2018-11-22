clear all ; close all ; 
EEG = pop_loadbv('C:\shared\raw_eeg\MONG_01_RB','MONG_01_RB_FIX_BOX.vhdr') ;
EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
rlabs = {EEG.chanlocs.labels} ; 

cd c:/shared/ute ; mongs=dir('*') ; mongs(1:2) = [] ; 
elecorder = {'FP1','FPZ','FP2','AF8','AF4','GND','AF3','AF7','F7','F5','F3','F1','FZ','F2','F4','F6','F8','FT10','FT8','FC6','FC4','FC2','REF','FC1','FC3','FC5','FT7','FT9',...
    'T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP10','TP8','CP6','CP4','CP2','CPZ','CP1','CP3','CP5','TP7','TP9','P7','P5','P3','P1','PZ','P2','P4','P6','P8',...
    'PO8','PO4','POZ','PO3','PO7','O1','OZ','O2'} ; 
% get the electrode indices
for i=1:length(rlabs) ; 
    if ~strcmp(rlabs{i},'ECG') ; 
        eeginds(i) = find(strcmpi(rlabs{i},elecorder)) ;    
    end
end

for mng=1:length(mongs) ;
    disp(mng) ; 
    cd(['c:/shared/ute/',mongs(mng).name,'/labs']) ; ls 
    allLabs = dir('*gz') ; labnames = {} ; 
    for i=1:length(allLabs) ;
        labsi = load_untouch_nii(allLabs(i).name) ;
        [lablocs(i,1),lablocs(i,2),lablocs(i,3)] = centmass3(labsi.img) ;
        labnames{i} = strrep(allLabs(i).name,'.label.nii.gz','') ; 
    end
    cd(['c:/shared/ute/',mongs(mng).name]) ; 
    coords = load('mricoords.mat') ; coords = coords.mricoords ; 
    %elabs = load('elabels.mat') ; elabs = elabs.elabels ; allelabs(mng,:) = elabs ; 
    %cskinter = load_untouch_nii('cskinter.nii.gz') ; cskinter = cskinter.img ; 
    %for i=1:63 ; [elocs(i,1),elocs(i,2),elocs(i,3)] = centmass3(cskinter==i) ; end
    
    for i=1:size(lablocs,1)
        for j=1:size(coords,2)
            distmat(mng,i,j) = sum(sqrt((lablocs(i,1)-coords(1,j)).^2 + (lablocs(i,2)-coords(2,j)).^2 + (lablocs(i,3)-coords(3,j)).^2)) ;       
        end
    end
end


mstdmat = squeeze(std(distmat,0,1)) ; 
mdistmat = squeeze(mean(distmat,1)) ; 
plot(reshape(mdistmat,[1,65*68]),reshape(mstdmat,[1,65*68]),'.','Color',[0,0,0]) ; title(['r = ',num2str(corr2(reshape(mstdmat,[1,65*68]),reshape(mdistmat,[1,65*68])))]) ; 
xlabel('distance (mm)') ; ylabel('std(mm)') ; 

for i=1:68 ; 
topod = squeeze(mdistmat(i,eeginds(find(eeginds~=0)))) ; 
fullvec = zeros(1,64) ; fullvec(1:31) = topod(1:31) ; fullvec(33:64) = topod(32:end) ; 
subplot(7,10,i) ; 
topoplot(fullvec,EEG.chanlocs,'maplimits',[0,30]) ; title(labnames{i}) ; colorbar
end

mdistmat = squeeze(mean(distmat,1)) ; clear fullvecs ; 
for s=1:8
for i=1:68 ; 
topod = squeeze(distmat(s,i,eeginds(find(eeginds~=0)))) ; 
fullvec = zeros(1,64) ; fullvec(1:31) = topod(1:31) ; fullvec(33:64) = topod(32:end) ; 
%subplot(7,10,i) ; 
%topoplot(fullvec,EEG.chanlocs,'maplimits',[0,100]) ; title(labnames{i}) ;
fullvecs(s,i,:) = fullvec ; 
end
end
for i=1:8 ; subplot(2,4,i) ; topoplot(squeeze(fullvecs(i,22,:)),EEG.chanlocs,'maplimits',[0,60]) ; title(['subject ',num2str(i)]) ; end

imagesc(squeeze(mean(distmat,1)))
set(gca,'YTick',1:68,'YTickLabel',labnames) ; 

fhandle = figure('Position',[10,-10,1600,900]) ; 
imagesc(mdistmat) ; hline(0.5:1:68,'k') ; vline(0.5:1:65,'k') ;  
set(gca,'YTick',1:68,'YTickLabel',labnames,'XTick',[]) ; axis xy

fhandle = figure('Position',[10,-10,1600,900]) ; 
plot(1:65) ; 
set(gca,'XTick',1:65,'XTickLabel',elecorder) ; rotateticklabel(gca,45) ; 

%{
stddists = (squeeze(std(distmat,[],1))) ;

for i=60 ; 
    figure,
    [sv,si] = sort(stddists(:,i)) ; subplot(2,1,1) ;
    bar(squeeze(stddists(si,i))) ; 
    set(gca,'XTick',1:68,'XTickLabel',labnames(si)) ;
    rotateticklabel(gca,45) ; ylabel('std(mm)') ; 
    suptitle(['standard deviation of distance between electrode ',baselabs{i},' and all cortical regions']) ; 
end

imagesc(stddists) ; 
set(gca,'YTick',1:68,'YTickLabel',labnames,'XTick',[]) ; axis xy
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') - [3 0,0])

plot(1:63) ; set(gca,'XTick',1:63,'XTickLabel',baselabs) ; 
rotateticklabel(gca,45)
%}










