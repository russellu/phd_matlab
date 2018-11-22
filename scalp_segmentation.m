%cd c:/shared/t1_ute ; clear all ; close all ; 
%subs = {'felix.nii.gz','MONG_01_RB.nii.gz','MONG_02_DP.nii.gz','MONG_03_CG.nii.gz','MONG_04_AB.nii.gz','MONG_05_SG.nii.gz','MONG_06_TS.nii.gz','nic.nii.gz'} ; 
ls 
padamt = 80 ; 
sphsize = 150 ; 
vecl = 148 ;
cd c:/shared/tests2
subs = {'t1_in_ute.nii.gz'} ; 

for sbb=1:length(subs)
disp(subs(sbb)) ; 
t1 = load_untouch_nii(subs{sbb}) ; 
t1img = t1.img ; 
flat1 = reshape(t1img,[size(t1img,1),size(t1img,2)*size(t1img,3)]) ; 
[kc,km] = kmeans(round(mat2gray(flat1)*255),10) ; 
kimg = reshape(km,size(t1img)) ; 

binimg = kimg>1 ; bwbin = bwconncomp(binimg) ; px = bwbin.PixelIdxList ; 
binzeros = zeros(size(binimg)) ; binzeros(px{1}) = 1 ; binimg = binzeros ; 

padbin = pad3d(binimg,padamt) ; 
prevdil = padbin ; 
for i=1:8 ; 
    dilmaski = imdilate(prevdil,strel(ones(3,3,3))) ; 
    clayer = (dilmaski - prevdil) > 0 ; 
    intensitylayers(clayer==1) = i ; 
    if i==1
        prevdil = dilmaski ;
    else
        prevdil = (prevdil + dilmaski) > 0 ; 
    end
end
padbin = prevdil ; 
[cx,cy,cz] = centmass3(padbin) ; 
xsi = size(padbin,1) ; ysi = size(padbin,2) ; zsi = size(padbin,3) ; 
[xg,yg,zg] = ndgrid(-cx:xsi-cx-1,-cy:ysi-cy-1,-cz:zsi-cz-1) ; 
sphere = sqrt(xg.^2 + yg.^2 + zg.^2) < sphsize ; sphere = imdilate(sphere,strel(ones(3,3,3)))-sphere ; 
sphereinds = find(sphere==1) ; 
[sphx,sphy,sphz] = ind2sub(size(padbin),sphereinds) ; % indices of the shell
xdiffs = sphx-cx ; ydiffs = sphy-cy ; zdiffs = sphz-cz ; 
% draw the vectors from the center of mass to the bounding points on the sphere:
diffvecs = [xdiffs,ydiffs,zdiffs]' ; 
unitdiffs = diffvecs./repmat(sqrt(sum(diffvecs.^2,1)),[3,1]) ; % normalized difference vectors (direction vectors)
unitreps = repmat(unitdiffs,[1,1,vecl]) ;
lines = zeros(1,1,vecl) ; lines(1,1,:) = 1:vecl ; lines = repmat(lines,[3,size(unitreps,2),1]) ; 
multlines = lines.*unitreps ; clear multinds
multinds(1,:,:) = floor(multlines(1,:,:)+cx) ; multinds(2,:,:) = floor(multlines(2,:,:)+cy) ; multinds(3,:,:) = floor(multlines(3,:,:)+cz) ; 
% multiply the lines by the brain, and get the highest index?
clear lineinds 
for i=1:size(multinds,3)
   lineinds(i,:) = sub2ind(size(padbin),squeeze(multinds(1,:,i)),squeeze(multinds(2,:,i)),squeeze(multinds(3,:,i))) ;
end
clearbrain = zeros(size(padbin)) ; 
clearbrain(lineinds(:,:)) = 1 ; 
%%% make the rbrain:
bininds = find(padbin==1) ; 
[bx,by,bz] = ind2sub(size(padbin),bininds) ; dbx = bx-cx ; dby = by-cy ; dbz = bz-cz ; 
[theta,phi,rho] = cart2sph(dbx,dby,dbz) ; 
rhobrain = zeros(size(padbin)) ; 
rhobrain(bininds) = rho ; 
clear rholines
for i=1:size(lineinds,1)
   rholines(i,:) = rhobrain(lineinds(i,:)) ;  
end
[sv,si] = sort(rholines,1,'descend') ; 
maxinds = si(1:20,:) ; 
clear scalpinds ; 
for i=1:size(lineinds,2) ; scalpinds(:,i) = lineinds(maxinds(:,i),i) ; end
%maxlines = max(rholines,[],1) ; 
%clear maxinds
%for i=1:size(rholines,2) ; maxinds(i) = lineinds(find(rholines(:,i)==maxlines(i),1),i) ; end

scalpbrain = zeros(size(padbin)) ; 
scalpbrain(scalpinds(:)) = 1 ; 
scalpres = scalpbrain(padamt:size(scalpbrain,1)-padamt-1,padamt:size(scalpbrain,2)-padamt-1,padamt:size(scalpbrain,3)-padamt-1) ; 
%filtimg = medfilt3(imdilate(scalpres,strel(ones(3,3,3)))) ; 

t1.img = scalpres.*(binimg==0) ; filtimg = t1.img ; 
trimmed = zeros(size(filtimg)) ; 
for i=1:size(filtimg,3) ;
   bw = bwconncomp(squeeze(filtimg(:,:,i))) ;
   if ~isempty(bw.PixelIdxList)
       inds = bw.PixelIdxList{1} ; 
       img = zeros(size(trimmed,1),size(trimmed,2)) ; 
       img(inds) = 1 ; 
       trimmed(:,:,i) = img ; 
   end
end
%t1.img = trimmed ; save_untouch_nii(t1,'trimmed.nii.gz') ; 

%trimlayers = (imdilate(trimmed,strel(ones(3,3,3)))-trimmed) ; 
%t1.img = trimlayers ; save_untouch_nii(t1,'trimmed.nii.gz') ; 
%bcomps = bwconncomp(trimlayers) ; 
%trimmed = medfilt3(trimmed) ; 
trimmed=  imdilate(trimmed,strel(ones(3,3,3))) ; 

[cx,cy,cz] = centmass3(trimmed) ; 
centind = sub2ind(size(trimmed),cx,cy,cz) ; 
conncomps = bwconncomp(trimmed==0) ; 
pixlist = conncomps.PixelIdxList ; 
% GET THE COMPONENT that contains the center of mass
for i=1:length(pixlist)
    if ismember(centind,pixlist{i}) ;
        braininds = pixlist{i} ; 
        break
    end
    disp(i) ; 
end
final = zeros(size(trimmed)) ; 
final(braininds) = 1 ; t1.img = imdilate(final,strel(ones(3,3,3))) ; 
save_untouch_nii(t1,['scalp_',subs{sbb}]) ;

end
