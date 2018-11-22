clear all ; close all ; 
cd c:/shared/utef/ ;
mongs=dir('*') ; mongs(1:2) = [] ; 
for m=3%:length(mongs) ; 
cd(['C:\shared\meanelecs\',mongs(m).name]) ; ls ; 
regs = dir('real_*') ; 
for reg = 1:length(regs)
    nii = load_untouch_nii(regs(reg).name) ; 
    x = medfilt3(nii.img) ; 
    zs = find(x==0) ; 
    x = x - min(min(min(x))) ;  
    x(zs) = 0 ; 
    diffx = abs(x-imfilter(x,fspecial('gaussian',3,3))) ;
    resx = uint8(mat2gray(reshape(diffx,[1,numel(x)]))*255) ; 
    [kc,km] = kmeanscustom(resx,2) ; 
    resk = reshape(km,size(x)) ; 
    x(resk==2) = 0 ; 
    allniis(reg,:,:,:) = x ; 
end
for i=1:65 ; 
    %figure ;
    x = squeeze(allniis(i,:,:,:)) ; 
    xres = uint8(mat2gray(reshape(x,[1,numel(x)]))*255) ; 
    [kc,km] = kmeanscustom(xres,5) ;
    x = reshape(km,size(x)) ; 
    newx(i,:,:,:) = x==5 ; 
    %for j=1:21 ; 
    %    subplot(4,6,j) ;      
    %    imagesc(squeeze(x(:,:,j))) ;
    %end  
end

cd(['C:\shared\utef\',mongs(m).name]) ; ls ; 
regute = load_untouch_nii('res_ute.nii.gz') ; 
zrute = zeros(size(regute.img)) ;
esize = 10 ; 
coords = load('mricoords') ; coords = coords.mricoords ; coords = coords + esize ; 
padz = pad3d(zrute,esize) ; 

for reg = 1:size(newx,1)
   padz(coords(1,reg)-esize:coords(1,reg)+esize,coords(2,reg)-esize:coords(2,reg)+esize,coords(3,reg)-esize:coords(3,reg)+esize) = squeeze(newx(reg,:,:,:)) ; 
end
regute.img = padz(esize:end-(esize+1),esize:end-(esize+1),esize:end-(esize+1)) ; 
save_untouch_nii(regute,'regute.nii.gz') ; 


end
