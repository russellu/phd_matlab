clear all  ;close all  ;
cd c:/shared/ute/MONG_01_RB ; ls ; 

ute = load_untouch_nii('t1_in_ute.nii.gz') ; 
rpial = load_untouch_nii('right_pial_in_ute.nii.gz') ; 
lpial = load_untouch_nii('left_pial_in_ute.nii.gz') ; 

padamt = 80 ; 

bothpial = rpial.img + lpial.img ; clear padbothpial
for i=1:size(bothpial,4) ; padbothpial(:,:,:,i) = pad3d(bothpial(:,:,:,i),padamt) ; end
pialmask = sum(abs(padbothpial),4) > 0 ; 
stdmask = load_untouch_nii('stdmask_in_ute.nii.gz') ; maskimg = pad3d(stdmask.img,padamt) ; 
pialmask = pialmask.*maskimg ; 
locs = load_untouch_nii('color_locs.nii.gz') ; origlocs = locs.img ; %origlocs = imdilate(origlocs,strel(ones(3,3,3))) ;
elecorder = load('elecorder.mat') ; elecorder = elecorder.elecorder ; 
locimg = pad3d(locs.img,padamt) ; % pad the locations image

projbrain = zeros(size(pialmask)) ; 
cortinds = find(pialmask==1) ; 
[cortx,corty,cortz] = ind2sub(size(pialmask),cortinds) ; 
for i=1:length(cortx)
    normalv = squeeze(padbothpial(cortx(i),corty(i),cortz(i),:)) ; 
    projvec = -90:90 ; projx = round(projvec*normalv(1)) ; projy = round(projvec*normalv(2)) ; projz = round(projvec*normalv(3)) ; 
    pvec = projvec ; pvec(projvec<0) = 1 ; pvec(projvec>=0) = -1 ; 
    dists = (90-sqrt(projx.^2 + projy.^2 + projz.^2)).*pvec ; 
    projx = projx + cortx(i) ; projy = projy + corty(i) ; projz = projz + cortz(i) ; 
    inds = sub2ind(size(projbrain),projx,projy,projz) ; 
    projbrain(inds) = projbrain(inds) + dists ;
end

ute.img = projbrain(padamt:size(projbrain,1)-padamt-1,padamt:size(projbrain,2)-padamt-1,padamt:size(projbrain,3)-padamt-1) ; 
save_untouch_nii(ute,'projbrain.nii.gz') ; 

for i=1:65 ; 
    locsi = find(origlocs==i) ;  
    intensities(i) = sum(ute.img(locsi)) ; 
end
figure,bar(abs(intensities)) ; set(gca,'XTick',1:65,'XTickLabel',elecorder) ; 
postelecs = elecorder(49:end) ; 
clear indies
for i=1:length(postelecs)  
    indies(i) = find(strcmpi(postelecs{i},elabs)) ; 
end

mersp = squeeze(mean(mean(mean(ersp,1),2),4)) ; 
modul = abs(squeeze(mean(mean(mersp(indies,freqs>10 & freqs<25,times>0 & times<10),2),3))) ; 
figure,plot(modul,abs(intensities(49:end)),'o') ;  title(num2str(corr2(modul,intensities(49:end)'))) ; 
