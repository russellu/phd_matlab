clear all ; close all ; 

cd c:/shared/utef/ ;
mongs=dir('*') ; mongs(1:2) = [] ; 

for m=1:length(mongs) ; 
    
cd(['C:\shared\utef\',mongs(m).name]) ; ls ; 

%locs = load_untouch_nii('locs.nii.gz') ; 
rute = load_untouch_nii('res_ute.nii.gz') ; ruteorig = double(rute.img) ; 
ruteorig = ruteorig-imfilter(ruteorig,fspecial('gaussian',80,40)) ; 
mask = load_untouch_nii('ants_ute/warpmask.nii.gz') ; mask.img = medfilt3(double(mask.img)) ; 
outermask = load_untouch_nii('ants_ute/warpmask_outer.nii.gz') ; dilouter = imdilate(outermask.img,strel(ones(21,21,21))) ; 
maskimg = mask.img>0 ; dilmask = imdilate(maskimg,strel(ones(21,21,21))) ; 
shell = dilmask-maskimg ; %shell = imdilate(shell,strel(ones(9,9,9))) ; % shell = imerode(shell,strel(ones(3,3,3))) ;  
shellrute = ruteorig.*shell.*dilouter ; 
rute.img = shellrute ; save_untouch_nii(rute,'shellrute.nii.gz')  ;

shell = dilmask-maskimg ; shellrute = ruteorig.*shell ; 
esize = 10 ; 
padmask = imerode(pad3d(maskimg,esize),strel(ones(17,17,17))) ; 
shellrute = pad3d(shellrute,esize) ; 
%shellrute = pad3d(ruteorig,esize) ; 
%shellrute(padmask==1) = mean(mean(mean(shellrute))) ;

coords = load('mricoords') ; coords = coords.mricoords ; coords = coords + esize ; 
clear boxes ; 
for i=1:size(coords,2)
    boxes(i,:,:,:) = medfilt3(shellrute(coords(1,i)-esize:coords(1,i)+esize,coords(2,i)-esize:coords(2,i)+esize,coords(3,i)-esize:coords(3,i)+esize)) ;   
end
resbox = reshape(boxes,[1,numel(boxes)])  ;
[kc,km] = kmeans(uint8(mat2gray(resbox)*255),3) ; 
%boxes = reshape(km,size(boxes)) ; 
allboxes(m,:,:,:,:) = boxes ; 
end
%for i=1:65 ; figure ; for j=1:21 ; subplot(3,7,j) ; imagesc(squeeze(boxes(i,:,:,j)),[0,550]) ; end ; end
%for i=1:65 ; figure ; for j=1:21 ; subplot(3,7,j) ; imagesc(squeeze(mboxes(i,:,:,j))) ; end ; end
%ab = (squeeze(allboxes(1,1,:,:,:))) ; 
%[kc,km] = kmeans(uint8(mat2gray(ab(ab~=0))*255),2) ;
%for i=1:65 ; figure ; for j=1:19 ; subplot(4,5,j) ; imagesc(squeeze(allboxes(1,i,:,:,j))) ; end ; end 

resbox = reshape(boxes,[1,numel(boxes)])  ;
[kc,km] = kmeans(uint8(mat2gray(resbox)*255),3) ; 
boxes = reshape(km,size(boxes)) ; 

cd c:/shared/ ; mkdir meanelecs 
cd meanelecs ; 
for s=1:length(mongs)
    mkdir(mongs(s).name) ; cd(mongs(s).name) ; 
    for i=1:65 ;
        if i>=10 ; stri = num2str(i) ; else stri = ['0',num2str(i)] ; end
        m3 = squeeze(allboxes(s,i,:,:,:)) ; 
        save_nii(make_nii(m3),['real_elec_',stri,'.nii.gz']) ; 
        %save_nii(make_nii(single(squeeze(m3~=0))),['mask_elec_',stri,'.nii.gz']) ; 
    end
    cd ..
end





%{
cd c:/shared/ATLASES ; 
mnimask = load_untouch_nii('mnimask_final.nii.gz') ; 
outer_mask = mnimask.img - imerode(mnimask.img,strel(ones(3,3,3))) ; 
mnimask.img = outer_mask ; 
save_untouch_nii(mnimask,'mnimask_outer.nii.gz') ; 
%}


cd c:/shared/utef/cloud ; 
regute = load_untouch_nii('res_ute.nii.gz') ; 
zrute = zeros(size(regute.img)) ;
coords = load('mricoords') ; coords = coords.mricoords ; coords = coords + 20 ; 
padz = pad3d(zrute,20) ; 

cd c:/shared/elecregs ; 
regs = dir('regmean_3*') ; 
for reg = 1:length(regs)
    nii = load_untouch_nii(regs(reg).name) ; 
  figure
   for j=1:size(nii.img,3) ; subplot(5,6,j) ; imagesc(squeeze(nii.img(:,:,j)),[-170,100]) ; end   
   padz(coords(1,reg)-8:coords(1,reg)+8,coords(2,reg)-8:coords(2,reg)+8,coords(3,reg)-8:coords(3,reg)+8) = squeeze(nii.img) ; 
   allniis(reg,:,:,:) = nii.img ; 
end
regute.img = padz(20:end-21,20:end-21,20:end-21) ; 
cd c:/shared/utef/cloud ; save_untouch_nii(regute,'regute.nii.gz') ; 


for e=1:65 ; 
    figure,
    img = squeeze(allniis(e,:,:,:)) ; 
    clear centers radii
    for i=1:size(img,3)    
        [centers{i},radii{i}] = imfindcircles(squeeze(img(:,i,:)),[5,15],'ObjectPolarity','bright') ; 
        subplot(4,5,i) ; imagesc(squeeze(img(:,i,:))) ; 
        if ~isempty(centers{i}) ;
            hline(round(centers{i}(2))) ; vline(round(centers{i}(1))) ; 
        end  
    end
end








