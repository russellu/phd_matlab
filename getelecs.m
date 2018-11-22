clear all ; close all ; 
cd c:/shared/lastute/ ;
mongs=dir('*') ; mongs(1:2) = [] ; 

for lab=1:15
for m=6%:length(mongs) ; 
    cd(['c:\shared\lastute\',mongs(m).name]) ; ls ; 
    mask = load_untouch_nii('finalmask.nii.gz') ; 
    maskim = imdilate(medfilt3(mask.img>0),strel(ones(3,3,3))) ; 
    resute = load_untouch_nii('res_ute.nii.gz') ; 
    resim = resute.img ; 
    coords = load(['mricoords_',num2str(lab),'.mat']) ; coords = coords.mricoords ; 
    
    allcoords(lab,:,:) = coords  ;
    esize = 10 ; 
    padmask = pad3d(maskim,esize) ; padres = pad3d(resim,esize) ; coords = coords + esize ; 
    padres = padres - imfilter(padres,fspecial('gaussian',80,50)) ; 
    shell = imdilate(padmask,strel(ones(7,7,7)))-padmask ; 
    shellrute = shell.*padres ; 
    zpad = zeros(size(shellrute)) ; zpad2 = zeros(size(shellrute)) ; 
    for i=1:size(coords,2)
        x = shellrute(coords(1,i)-esize:coords(1,i)+esize,coords(2,i)-esize:coords(2,i)+esize,coords(3,i)-esize:coords(3,i)+esize) ; 
        inds = find(x==0) ; 
        x = x-min(min(min(x))) ;     
        x(inds) = 0 ; 
        kinds = find(x~=0) ; kvals = uint8(mat2gray(x(kinds))*255) ; 
        [kc,km] = kmeanscustom(kvals,5) ; 
        x(kinds) = km ; 
        elecs(i,:,:,:) = x ; 
        zpad(coords(1,i)-esize:coords(1,i)+esize,coords(2,i)-esize:coords(2,i)+esize,coords(3,i)-esize:coords(3,i)+esize) = (x>4).*i ; 
    end
    resute.img = zpad(esize:end-(esize+1),esize:end-(esize+1),esize:end-(esize+1)) ; 
    allzpads(lab,:,:,:) = resute.img ;     
end
end
for i=1:65 ; figure ; for j=1:11 ; subplot(4,6,j) ; imagesc(squeeze(elecs(i,:,:,j)),[0,5]) ; end ; end

clear cx cy cz
for i=1:65
    for j=1:15
        [cx(i,j),cy(i,j),cz(i,j)] = centmass3(squeeze(allzpads(j,:,:,:))==i) ;         
    end
end

clear segs
segs(:,1,:) = cx' ; segs(:,2,:) = cy' ; segs(:,3,:) = cz' ; 
labstd = squeeze(std(allcoords,0,1)) ; 
segstd = squeeze(std(segs,0,1)) ; 
elecorder = {'FP1','FPZ','FP2','AF8','AF4','GND','AF3','AF7','F7','F5','F3','F1','FZ','F2','F4','F6','F8','FT10','FT8','FC6','FC4','FC2','REF','FC1','FC3','FC5','FT7','FT9',...
    'T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP10','TP8','CP6','CP4','CP2','CPZ','CP1','CP3','CP5','TP7','TP9','P7','P5','P3','P1','PZ','P2','P4','P6','P8',...
    'PO8','PO4','POZ','PO3','PO7','O1','OZ','O2'} ; 
bads = {'FT9','FT10','TP9','TP10'} ; badinds = [18,28,38,48] ; 
labstd(:,badinds) = [] ; 

segstd(:,badinds) = [] ; 
subplot(1,2,1) ; bar(sort(mean(segstd))) ; title(num2str(mean(mean(segstd)))) ; ylim([0,3]) ; 
subplot(1,2,2) ; bar(sort(mean(labstd))) ; title(num2str(mean(mean(labstd)))) ; ylim([0,3]) ; 
