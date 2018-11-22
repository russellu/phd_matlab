clear all ; close all ; 
cd c:/shared/lastute ; 
subs=dir('*') ; 
subs(1:2) = [] ; 
for s=1:length(subs)
cd(['c:/shared/lastute/',subs(s).name]) ;  
coords = dir('mricoords_*') ; 
clear allmricoords
for coord=1:length(coords)
    mricoords = load(coords(coord).name) ; 
    mricoords = mricoords.mricoords ; 
    allmricoords(coord,:,:) = mricoords ; 
end

%elecorder = load('elecorder.mat') ; elecorder = elecorder.elecorder ; 
allmricoords(:,:,[18,28,38,48]) = [] ; 
mask = load_untouch_nii('finalmask.nii.gz') ; 
ute = load_untouch_nii('res_ute.nii.gz') ; 
padute = pad3d(ute.img,25) ; % padute = padute-imfilter(padute,fspecial('gaussian',80,40)) ; 
padmask = pad3d(mask.img,25) ; 
allmricoords = allmricoords + 25 ; 
shell = imdilate(padmask,strel(ones(5,5,5))) - padmask ; 
shellute = shell.*padute ; 
esize = 10 ; clear eboxes ; 
for i=1:size(allmricoords,1)
    for j=1:size(allmricoords,3)  
        eboxes(i,j,:,:,:) = shellute(allmricoords(i,1,j)-esize:allmricoords(i,1,j)+esize,...
                                     allmricoords(i,2,j)-esize:allmricoords(i,2,j)+esize,...
                                     allmricoords(i,3,j)-esize:allmricoords(i,3,j)+esize) ; 
    end
end

newbrains = zeros(15,size(padute,1),size(padute,2),size(padute,3)) ; 
squareboxes = zeros(size(padute)) ; 
for i=1:size(eboxes,1) ; disp(i)
    for j=1:size(eboxes,2) ; disp(j) ; 
        m3 = medfilt3(squeeze(eboxes(i,j,:,:,:))) ; 
        inds = find(m3 ~= 0) ;
        vals = m3(inds) ; 
        [~,km] = kmeanscustom(uint8(mat2gray(vals)*255),4) ; 
        m3(inds) = km ; 
        newbrains(i,allmricoords(i,1,j)-esize:allmricoords(i,1,j)+esize,...
                    allmricoords(i,2,j)-esize:allmricoords(i,2,j)+esize,...
                    allmricoords(i,3,j)-esize:allmricoords(i,3,j)+esize) = (m3==4).*j ;  
        if i==1
            squareboxes(allmricoords(i,1,j)-esize:allmricoords(i,1,j)+esize,...
                    allmricoords(i,2,j)-esize:allmricoords(i,2,j)+esize,...
                    allmricoords(i,3,j)-esize:allmricoords(i,3,j)+esize) = 1 ;  
        end
    end
end    

resnew = newbrains(:,25:end-26,25:end-26,25:end-26) ; 
for i=1:15 ; disp(['centmass: ',num2str(i)]) ; 
    for j=1:size(eboxes,2)
        [cx,cy,cz] = centmass3(squeeze(resnew(i,:,:,:))==j) ; 
        allcms(s,i,j,:) = [cx,cy,cz] ; 
    end
end
end

mcms = squeeze(mean(allcms,2)) ; 
for s=1:length(subs)
cd(['c:/shared/lastute/',subs(s).name]) ;  
smcm = squeeze(mcms(s,:,:)) ; 
save('meanlocs','smcm') ; 
end













