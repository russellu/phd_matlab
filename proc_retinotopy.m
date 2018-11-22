cd c:/shared/retinotopic ; ls ; clear all ; close all ; 
% process retinotopic data
stimvecs = dir('stimvec*') ; 
for stimv = 1 % get the stimulus log file (all 3 are the same, however)
    svi = load(stimvecs(stimv).name) ; 
    svi = svi.stimvec ;  
    xys(1,:) = squeeze(svi(:,1)) ; 
    xys(2,:) = squeeze(svi(:,2)) ; 
    t = squeeze(svi(:,3)) ; 
end
xys(:,1) = [] ; t(1) = [] ;
TR = 0.41 ; 
trs = round(t/TR) ;
raws = dir('reg*') ; 
for r=1:3
    nii = load_untouch_nii(raws(r).name) ; 
    nimg = nii.img ; 
    % epoch
    baselength = floor(2/TR) ; tasklength = floor(8/TR) ; 
    for i=1:length(t)
        epochs(:,:,:,i,:) = nimg(:,:,:,trs(i)-baselength:trs(i)+tasklength) ;     
    end

    diffepochs = (squeeze(mean(epochs(:,:,:,:,13:end),5))-squeeze(mean(epochs(:,:,:,:,1:baselength),5)))./squeeze(mean(epochs(:,:,:,:,1:baselength),5)) ; 
    alldiffepochs(r,:,:,:,:) = diffepochs ; 

end
%diffepochs = squeeze(mean(alldiffepochs)) ; 
diffepochs = double((squeeze(mean(epochs(:,:,:,:,13:end),5))-squeeze(mean(epochs(:,:,:,:,1:baselength),5)))./squeeze(mean(epochs(:,:,:,:,1:baselength),5))) ; 
rads = sqrt((xys(1,:)-mean(xys(1,:))).^2 + (xys(2,:)-mean(xys(2,:))).^2) ; 
goodrads = find(rads<300) ; 
diffepochs = diffepochs(:,:,:,goodrads) ; 
xbrain = rads(goodrads) ; 
cvol = voxcorr(diffepochs,xbrain) ; 
figure,for i=50:64 ; subplot(3,5,i-49) ; imagesc(rot90(squeeze(cvol(:,i,:))),[-.3,.3]); title(i) ; end 
figure,for i=1:18 ; subplot(4,5,i) ; imagesc((squeeze(cvol(:,:,i))),[-.3,.3]); title(i) ; end 
inds = find(cvol>-.2) ; [cx,cy,cz] = ind2sub(size(cvol),inds) ; 
[xq,yq] = meshgrid(min(xys(1,:)):max(xys(1,:)),min(xys(2,:)):max(xys(2,:))) ; 
for i=1:length(cx) ; 
    rets  = griddata(xys(1,goodrads),xys(2,goodrads),squeeze(diffepochs(cx(i),cy(i),cz(i),:)),xq,yq,'v4') ;     
    figure ; imagesc(rets) ; 
    meanrets(i,:,:) = rets ; 
end
figure,imagesc(squeeze(mean(meanrets))) ; 
f1 = double(squeeze(nimg(:,:,:,1))>200) ; 
for i=1:18
subplot(3,6,i),plotoverlayIntensity2D(squeeze(nimg(:,:,i,1).*f1(:,:,i)),sqrt(abs(cvol(:,:,i)).*f1(:,:,i)),mat2gray(cvol(:,:,i).*f1(:,:,i)),270)
end




