clear all ; close all ; 
cd c:/shared/lastute ; 
subs=dir('*') ; 
subs(1:2) = [] ;     
for s=1:length(subs)
   cd(['c:/shared/lastute/',subs(s).name]) ;  
    t1 = load_untouch_nii('t1_in_ute.nii.gz') ; 
    EEG = pop_loadset('C:\shared\badger_eeg\alex\1hz_preproc_retino_allstims_01_Pulse Artifact Correction.set') ; 
    elabs = {EEG.chanlocs.labels} ;
    elecorder = {'FP1','FPZ','FP2','AF8','AF4','GND','AF3','AF7','F7','F5','F3','F1','FZ','F2','F4','F6','F8','FT10','FT8','FC6','FC4','FC2','REF','FC1','FC3','FC5','FT7','FT9',...
        'T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP10','TP8','CP6','CP4','CP2','CPZ','CP1','CP3','CP5','TP7','TP9','P7','P5','P3','P1','PZ','P2','P4','P6','P8',...
        'PO8','PO4','POZ','PO3','PO7','O1','OZ','O2'} ; 

    coords = load('meanlocs') ; coords = coords.smcm ; 
    cd fast ; 
    gm = load_untouch_nii('fast_seg_1.nii.gz')  ;
    csf = load_untouch_nii('fast_seg_0.nii.gz') ; 
    wm = load_untouch_nii('fast_seg_2.nii.gz') ; 
    xs = size(t1.img,1) ; ys = size(t1.img,2) ; zs = size(t1.img,3) ;
    for i=1:length(coords) ; disp(i) ; 
        cxyz = coords(i,:) ;     
        [gx,gy,gz] = ndgrid(-cxyz(1):xs-cxyz(1)-1,-cxyz(2):ys-cxyz(2)-1,-cxyz(3):zs-cxyz(3)-1) ; 
        for r=2:40 ;
            sphere = (sqrt(gx.^2 + gy.^2 + gz.^2) < r) & (sqrt(gx.^2 + gy.^2 + gz.^2) > r-1) ;  
            sumgm(s,i,r) = sum(sum(sum(double(sphere).*double(gm.img)))) ; 
            sumcsf(s,i,r) = sum(sum(sum(double(sphere).*double(csf.img)))) ; 
            sumwm(s,i,r) = sum(sum(sum(double(sphere).*double(wm.img)))) ; 
        end
    end
    disp(subs(s).name) ; 
end


o2 = squeeze(sumcsf(:,61,:)) ; 

newsphere = sqrt(gx.^2 + gy.^2 + gz.^2).*(sqrt(gx.^2 + gy.^2 + gz.^2)<35) ; 



plot(squeeze(mean(mean(sumgm,1),2)),'k','LineWidth',2) ; hold on ; plot(squeeze(mean(mean(sumcsf,1),2)),'b','LineWidth',2) ;plot(squeeze(mean(mean(sumwm,1),2)),'r') ;
xlabel('mm') ; ylabel('# voxels') ; legend({'gray matter','csf','white matter'}) ;

bar(squeeze(mean(mean(sumgm(:,:,20:end),1),3)));
    elecorder = {'FP1','FPZ','FP2','AF8','AF4','GND','AF3','AF7','F7','F5','F3','F1','FZ','F2','F4','F6','F8','FT10','FT8','FC6','FC4','FC2','REF','FC1','FC3','FC5','FT7','FT9',...
        'T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP10','TP8','CP6','CP4','CP2','CPZ','CP1','CP3','CP5','TP7','TP9','P7','P5','P3','P1','PZ','P2','P4','P6','P8',...
        'PO8','PO4','POZ','PO3','PO7','O1','OZ','O2'} ; 
    inds = 1:65 ; inds([18,28,38,48]) = [] ; 
    neworder = elecorder(inds);  
imagesc(squeeze(mean(sumgm(:,:,:),1))); set(gca,'YTick',1:61,'YTickLabel',neworder) ; xlabel('dist(mm)') ; hline(0:60,'k') ; 
%{
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
clear topodiffs
for i=1:length(revinds)
    if revinds(i) ~= 0 
        topogms(i) = sums(revinds(i)) ;     
    end
end
%}
    EEG = pop_loadset('C:\shared\badger_eeg\alex\1hz_preproc_retino_allstims_01_Pulse Artifact Correction.set') ; 
