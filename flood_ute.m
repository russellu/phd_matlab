cd c:/shared/tests2 ; ls ; clear all ; close all ; 
rute = load_untouch_nii('res_3.nii.gz') ; 
[kc,km] = kmeans(uint8(mat2gray(reshape(rute.img,[1,numel(rute.img)]))*255),2) ; 
kimg = reshape(km,size(rute.img)) > 1 ; 
padamt = 60 ; 
padded = pad3d(kimg,padamt) ; 
padk = pad3d(kimg,padamt) ; 
[cx,cy,cz] = centmass3(abs(padded)) ; [xg,yg,zg] = ndgrid(-cx:size(padded,1)-cx-1,-cy:size(padded,2)-cy-1,-cz:size(padded,3)-cz-1) ; 
[th,phi,rho] = cart2sph(xg,yg,zg) ; 
rhomask = rho.*(padk) ; 
sphere = sqrt(xg.^2 + yg.^2 + zg.^2) < 120 ; rhomask = rhomask.*sphere ; 
shell = sphere-imerode(sphere,strel(ones(3,3,3))) ; shellinds = find(shell==1) ; 
[sx,sy,sz] = ind2sub(size(shell),shellinds) ;
diffx = sx-cx ; diffy = sy-cy ; diffz = sz-cz ; difflengths = sqrt(diffx.^2 + diffy.^2 + diffz.^2) ; 
udiffx = diffx./difflengths ; udiffy = diffy./difflengths ; udiffz = diffz./difflengths ;
lines = repmat((1:round(mean(difflengths)))',[1,length(difflengths)]) ; 
clear xlines ylines zlines
for i=1:size(lines,1) ; 
    xlines(i,:) = floor(cx - lines(i,:).*udiffx')+1 ; 
    ylines(i,:) = floor(cy - lines(i,:).*udiffy')+1 ; 
    zlines(i,:) = floor(cz - lines(i,:).*udiffz')+1 ; 
end
clear inds vals
for i=1:size(lines,1)
    inds(i,:) = sub2ind(size(shell),xlines(i,:),ylines(i,:),zlines(i,:)) ; 
    vals(i,:) = rhomask(inds(i,:)) ; 
end
[maxvals,maxinds] = max(vals,[],1) ; 
zbrain = zeros(size(padded)) ; 
for i=1:size(maxinds,2)
    indi = inds(maxinds(i),i) ;
    zbrain(indi) = 1 ; 
end
newbrain = zbrain(padamt:size(zbrain,1)-padamt-1,padamt:size(zbrain,2)-padamt-1,padamt:size(zbrain,3)-padamt-1) ; 
%%%% do a region growing approach (using ray tracing)
slicebrain = zeros(size(newbrain)) ; 
for i=1:size(newbrain,3) 
     slicei = squeeze(newbrain(:,:,i)) ; 
     for j=1:size(slicei,1) % for all rows
         filling = true ; fillcount = 1 ; 
         for k=1:size(slicei,2) % for all columns
             if slicei(j,k)==1 
                 fillcount = fillcount + 1 ;                 
             end
             if fillcount > 1 && slicei(j,k) == 0
                fillcount = 1 ; filling = ~filling ;  
             end
             if filling 
                 slicei(j,k) = 1 ;                  
             end                           
         end
     end   
     slicebrain(:,:,i) = slicei ; 
end
slicebrain2 = zeros(size(newbrain)) ; 
for i=1:size(newbrain,3) 
     slicei = squeeze(newbrain(:,:,i)) ; 
     for j=1:size(slicei,2) % for all rows
         filling = true ; fillcount = 1 ; 
         for k=1:size(slicei,1) % for all columns
             if slicei(k,j)==1 
                 fillcount = fillcount + 1 ;                 
             end
             if fillcount > 1 && slicei(j,k) == 0
                fillcount = 1 ; filling = ~filling ;  
             end
             if filling 
                 slicei(k,j) = 1 ;                  
             end                           
         end
     end   
     slicebrain2(:,:,i) = slicei ; 
end
slicebrain3 = zeros(size(newbrain)) ; 
for i=1:size(newbrain,3) 
     slicei = squeeze(newbrain(:,:,i)) ; 
     for j=size(slicei,2):-1:1 % for all rows
         filling = true ; fillcount = 1 ; 
         for k=size(slicei,1):-1:1 % for all columns
             if slicei(k,j)==1 
                 fillcount = fillcount + 1 ;                 
             end
             if fillcount > 1 && slicei(j,k) == 0
                fillcount = 1 ; filling = ~filling ;  
             end
             if filling 
                 slicei(k,j) = 1 ;                  
             end                           
         end
     end   
     slicebrain3(:,:,i) = slicei ; 
end
slicebrain4 = zeros(size(newbrain)) ; 
for i=1:size(newbrain,3) 
     slicei = squeeze(newbrain(:,:,i)) ; 
     for j=size(slicei,1):-1:1 % for all rows
         filling = true ; fillcount = 1 ; 
         for k=size(slicei,2):-1:1 % for all columns
             if slicei(j,k)==1 
                 fillcount = fillcount + 1 ;                 
             end
             if fillcount > 1 && slicei(j,k) == 0
                fillcount = 1 ; filling = ~filling ;  
             end
             if filling 
                 slicei(j,k) = 1 ;                  
             end                           
         end
     end   
     slicebrain4(:,:,i) = slicei ; 
end

sumslice = slicebrain+slicebrain2+slicebrain3+slicebrain4 ; 
binslice = medfilt3(sumslice)> 2 ; 
bw = bwconncomp(binslice) ; clear lengths ; 
pixlist = bw.PixelIdxList ; for i=1:length(pixlist) ; lengths(i) = length(pixlist{i}) ; end
compimg = zeros(size(slicebrain)) ; compimg(pixlist{find(lengths==max(lengths))}) = 1 ; 
compimg = ~compimg ; 
bw = bwconncomp(compimg) ; clear lengths ; 
pixlist = bw.PixelIdxList ; for i=1:length(pixlist) ; lengths(i) = length(pixlist{i}) ; end
compimg = zeros(size(slicebrain)) ; compimg(pixlist{find(lengths==max(lengths))}) = 1 ; 
mederode = medfilt3(imerode(compimg,strel(ones(3,3,3)))) ; 
rute.img = mederode ; save_untouch_nii(rute,'mederode.nii.gz') ; 


