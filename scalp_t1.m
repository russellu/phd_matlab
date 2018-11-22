%%%% get the projection vectors from the inside of the brain outwards and
%%%% run ICA on that...
cd c:/shared/tests2 ; ls ; clear all ; close all ; 
ute = load_untouch_nii('t1_in_ute.nii.gz') ; 
utim = ute.img ; padamt = 120 ; 
% resize the image
%utim2 = (utim(1:2:end,1:2:end,1:2:end) + utim(2:2:end,2:2:end,2:2:end))./2 ; 
padded = pad3d(utim,padamt) ; 
[cx,cy,cz] = centmass3(abs(padded)) ; [xg,yg,zg] = ndgrid(-cx:size(padded,1)-cx-1,-cy:size(padded,2)-cy-1,-cz:size(padded,3)-cz-1) ; 
[th,phi,rho] = cart2sph(xg,yg,zg) ; 
[kc,km] = kmeans(uint8(mat2gray(reshape(padded,[1,numel(padded)]))*255),5) ; 
kimg = (reshape(km,size(padded))) ; 
rhomask = rho.*(kimg>1) ; 
sphere = sqrt(xg.^2 + yg.^2 + zg.^2) < 190 ; 
shell = imdilate(sphere,strel(ones(3,3,3))) - sphere ; shellinds = find(shell==1) ; 
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
clear scalpinds ; scalpinds(1,:) = maxinds ; 
for i=2:8 
   scalpinds(i,:) = scalpinds(1,:) + i ;    
end
zbrain = zeros(size(shell)) ; 
for i=1:size(scalpinds,1)
    maxinds2 = scalpinds(i,:) ; 
    for j=1:length(maxinds2)
        indj = inds(maxinds2(j),j) ; 
        zbrain(indj) = i ;    
    end
end
newz = zbrain(padamt:size(zbrain,1)-padamt-1,padamt:size(zbrain,2)-padamt-1,padamt:size(zbrain,3)-padamt-1) ; 
ute.img = newz ; save_untouch_nii(ute,'intensity_layers.nii.gz') ; 

%{
newz = zbrain(padamt:size(zbrain,1)-padamt-1,padamt:size(zbrain,2)-padamt-1,padamt:size(zbrain,3)-padamt-1) ; 
[xg1,yg1,zg1] = ndgrid(1:size(newz,1),1:size(newz,2),1:size(newz,3)) ; 
[xg2,yg2,zg2] = ndgrid(1:.5:size(newz,1)+.5,1:.5:size(newz,2)+.5,1:.5:size(newz,3)+.5) ; 
P = [2 1 3];
xg1 = permute(xg1, P);
yg1 = permute(yg1, P);
zg1 = permute(zg1, P);
newz = permute(newz, P);
vq = interp3(xg1,yg1,zg1,newz,xg2,yg2,zg2,'nearest') ; 
%}
%{
for i=1:length(maxinds)
    indi = inds(maxinds(i),i) ; 
    zbrain(indi) = zbrain(indi) + 1 ;    
end
newbrain = zbrain(padamt:size(zbrain,1)-padamt-1,padamt:size(zbrain,2)-padamt-1,padamt:size(zbrain,3)-padamt-1) ; 
disp3d(newbrain); 
%}


