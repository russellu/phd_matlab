%%%% MATLAB script to process FMRI data from visual experiments
% gets the stimulus times from log files, extracts stimulus type and time
% at which it occured, and saves an ideal convolved with an HRF to the
% subject specific directory
clear all ; close all

stimnames{1} = 'unperturbed'  ;
stimnames{2} = 'contrast_5%'  ;
stimnames{3} = 'contrast_33%'  ;
stimnames{4} = 'plaid'  ;
stimnames{5} = 'rnd_10%'  ;
stimnames{6} = 'rnd_60%'  ;

subjects = {
 %   'alex'
 %   'charest'
    'esteban'
    'fabio'
    'gab'
    'gabriella'
    'genevieve'
    'gina'
  %  'guillaume'
    'jeremie'
    'julie'
 %   'katrine'
    'lisa'
    'marc'
    'marie'
   'mathieu'
    'maxime'
    'mingham'
    'patricia'
    'po'
    'russell'
    'sunachakan'
  %  'tah'
  %  'thititip'
    'vincent'
} ; 
MNI = load_nii('c:/shared/allfmris/sub_russell/ants/warped.nii.gz') ; mim = MNI.img ; 
a = load_nii('c:/shared/stats/0a.nii.gz') ; a= a.img ;mask= a>.3;  


for sub=1:size(subjects,1) ; 
    cd(['c:/shared/allfmris/sub_',subjects{sub},'/ants']) ;
    ls
    cur = load_nii('warpcorrs.nii.gz') ; 
    corrs(sub,:,:,:) = cur.img; 
   
    warps=dir('warp_t1*') ; 
    for i=1:size(warps,1)
       nii = load_nii(warps(i).name) ; 
       imgs(sub,i,:,:,:) = nii.img ; 
    end
   
    %a = load_nii('warped.nii.gz') ; 
    
    
    
end
mnis = imgs ; 
%{
mask = mim > mean(mean(mean(mim))) ; 
pmat = zeros(182,218,182) ; 
for i=1:182
    disp(i) ;
    for j=1:218
        for k=1:182
            if mask(i,j,k) ~= 0
            	pmat(i,j,k) = anova1(squeeze(imgs(:,[1,5,6],i,j,k)),[],'off') ;   
            end
        end
    end
end
  %}


contrast = load_nii('c:/shared/stats/invcontrast.nii.gz') ; stims(1,:,:,:) = contrast.img ;
rnd = load_nii('c:/shared/stats/rnd.nii.gz') ; stims(2,:,:,:) = rnd.img ;
plaid = load_nii('c:/shared/stats/invplaid.nii.gz') ; stims(3,:,:,:) = plaid.img ;

sthresh = stims > .95 ; 
mcorrs = squeeze(mean(corrs,1)) ; 
cthresh = mcorrs > .25 ; 
for i=1:3 ; rois(i,:,:,:) = squeeze(sthresh(i,:,:,:)).*cthresh ; end
for i=1:size(imgs,1) ; for j=1:3 ; for k=1:6 ; percs(i,j,k) = (sum(sum(sum(squeeze(imgs(i,k,:,:,:)).*squeeze(rois(j,:,:,:))))))./sum(sum(sum(squeeze(rois(j,:,:,:))))) ; end ; end ; end 
for i=1:3 ; subplot(2,2,i) ; errorbar(squeeze(mean(percs(:,i,[1,4]),1)),squeeze(std(percs(:,i,[1,4]),0,1))./sqrt(21)) ; end


for thresh=0.1:0.05:0.35
 vol = squeeze(sum(sum(sum(mcorrs > thresh)))) ; mask = mcorrs> thresh ; for i=1:size(imgs,1) ; for j=1:size(imgs,2) ; svals(i,j) = sum(sum(sum(squeeze(imgs(i,j,:,:,:)).*mask))) ; end ; end
figure,subplot(2,2,1) ; errorbar(mean(svals(:,[1,3,2]),1),std(svals(:,[1,3,2]),0,1)./sqrt(21)) ; subplot(2,2,2); errorbar(mean(svals(:,[1,5,6]),1),std(svals(:,[1,5,6]),0,1)./sqrt(21)) ; subplot(2,2,3); errorbar(mean(svals(:,[1,4]),1),std(svals(:,[1,4]),0,1)./sqrt(21)) ; 
end
 



for i=1:3 ; rois(i,:,:,:) = squeeze(sthresh(i,:,:,:)).*cthresh ; end
for i=1:size(imgs,1) ; for j=1:3 ; for k=1:6 ; percs(i,j,k) = (sum(sum(sum(squeeze(imgs(i,k,:,:,:)).*squeeze(rois(j,:,:,:))))))./sum(sum(sum(squeeze(rois(j,:,:,:))))) ; end ; end ; end 
for i=1:3 ; subplot(2,2,i) ; errorbar(squeeze(mean(percs(:,i,[1,4]),1)),squeeze(std(percs(:,i,[1,4]),0,1))./sqrt(21)) ; end

tcount = 1 ; clear cmat
for t=0.1:0.05:0.35
    t
   cmask = mcorrs > t ; sums = sum(sum(sum(cmask))) ; 
   for i=1:size(imgs,1) ; for j=1:size(imgs,2) ; cmat(tcount,i,j) = squeeze(sum(sum(sum(cmask.*squeeze(imgs(i,j,:,:,:))))))./sums ; end ; end ; tcount = tcount + 1 ; 
end


clear pmats1 pmats2 pmats3 ; 
for i=1:6 ; pmats1(i) = anova1(squeeze(cmat(i,:,[1,3,2])),[],'off') ; end ; for i=1:6 ; pmats2(i) = anova1(squeeze(cmat(i,:,[1,5,6])),[],'off') ; end ; for i=1:6 ; pmats3(i) = anova1(squeeze(cmat(i,:,[1,4])),[],'off') ; end


subplot(3,1,1) ; bar(pmats1) ; set(gca,'XTick',1:26) ; set(gca,'XTickLabel',t) ; ylabel('p value') ; xlabel('correlation threshold (r)') ; title('differences due to contrast') ;
subplot(3,1,2) ; bar(pmats2) ; set(gca,'XTick',1:26) ; set(gca,'XTickLabel',t) ; ylabel('p value') ; xlabel('correlation threshold (r)') ; title('differences due to randomizaton') ;
subplot(3,1,3) ; bar(pmats3) ; set(gca,'XTick',1:26) ; set(gca,'XTickLabel',t) ; ylabel('p value') ; xlabel('correlation threshold (r)') ; title('differences due to plaid') ;

itvl = 62; i=1:14 ; i=1
img = fliplr(rot90(mat2gray(squeeze(mean(a(:,:,itvl),3)).^3))) ; 
im2 = (squeeze(mean(mean(corrs(i,:,:,itvl),1),4))) ; 
I = fliplr(rot90(uint8((mat2gray(im2)*255)))) ;
rgb = ind2rgb(gray2ind(I,255),jet(255)) ;
imshow(img, 'InitialMag', 'fit') ;
hold on ; h = imshow(rgb) ; hold off ; 
set(h, 'AlphaData', fliplr((rot90(squeeze(mean(mean(corrs(i,:,:,itvl),4),1)))))>.05) ;

%%% plot contrast and randomization f-values 
cd c:/shared/stats ; 
contrast = load_nii('fcontrast.nii.gz') ; contrast = contrast.img ; contrast = contrast.*(mcorrs>.1) ; 
rnd = load_nii('frnd.nii.gz') ; rnd = rnd.img ; rnd = rnd.*(mcorrs>.1) ;
fbrain = rnd ;
icount = 1 ; x = 25:160 ; y = 20:200 ; itincr = 3 ;
for itvl = 58:2:65
    subplot(2,3,icount) ; icount = icount + 1 ;
    img = fliplr(rot90(mat2gray(squeeze(mean(a(x,y,itvl:itvl+itincr),3)).^3))) ; 
    im2 = (squeeze(mean(fbrain(x,y,itvl:itvl+itincr),3))) ; 
    I = fliplr(rot90(uint8((mat2gray(im2)*255)))) ;
    rgb = ind2rgb(gray2ind(I,255),hot(255)) ;
    imshow(img, 'InitialMag', 'fit') ;
    hold on ; h = imshow(rgb) ; hold off ; 
    set(h, 'AlphaData', fliplr(abs(rot90(squeeze(mean(fbrain(x,y,itvl:itvl+itincr),3)))))>3) ;
end

%%% plot both on an rgb map

hg = load_nii('c:/shared/stats/rbrainhg.nii.gz') ; hg = hg.img ; beta = load_nii('c:/shared/stats/rbrainbeta.nii.gz') ; beta = beta.img ;
contrmask = contrast ; rndmask = rnd ; 
clear rgbs
rgbs(:,:,:,1) = contrmask>2 ; rgbs(:,:,:,3) = 0 ; rgbs(:,:,:,2) = rndmask>2 ; 
icount = 1 ; x = 25:160 ; y = 20:200 ; itincr = 3 ;
for itvl = 60
    subplot(1,1,1) ; icount = icount + 1 ;
    img = fliplr(rot90(mat2gray(squeeze(mean(a(x,y,itvl:itvl+itincr),3)).^3))) ; 
   % im2 = (squeeze(mean(contrast(x,y,itvl:itvl+itincr),3))) ; 
   for i=1:3
    rgb(:,:,i) = fliplr(rot90((squeeze(mean(rgbs(x,y,itvl:itvl+itincr,i),3))))) ; 
   end
% I = fliplr(rot90(uint8((mat2gray(im2)*255)))) ;
%    rgb = ind2rgb(gray2ind(I,255),hot(255)) ;
    imshow(img, 'InitialMag', 'fit') ;
    hold on ; h = imshow(mat2gray(rgb)) ; hold off ; 
    
    set(h, 'AlphaData', mean(rgb,3)>.2 ) ;
end

img = rot90(squeeze(mean(mat2gray(a(90,:,:)),1))) ; imcopy =img ; imcopy(110:120,:) = 1 ;
hold on ; imshow(img, 'InitialMag', 'fit') ;
h = imshow(imcopy) ; hold off ; 
set(h, 'AlphaData', imcopy) ;

%%% get the error bars from the single voxel data
mcorrs = squeeze(mean(corrs,1)) ; maskcorrs = mcorrs >.25 ;
hgthresh = hg > .9 ; bthresh = beta < -.9 ; 
bcount = 1 ; gcount = 1 ; 
clear hgs betas
for i=1:size(maskcorrs,1)
    for j=1:size(maskcorrs,2)
        for k=1:size(maskcorrs,3)
            if hgthresh(i,j,k)==1 && maskcorrs(i,j,k)==1
                hgs(gcount,:,:) = squeeze(mnis(:,:,i,j,k)) ; gcount = gcount + 1 ;
            end
            if bthresh(i,j,k)==1 && maskcorrs(i,j,k)==1
                betas(bcount,:,:) = squeeze(mnis(:,:,i,j,k)) ; bcount = bcount + 1 ;
            end
        end
    end
end
stims = [2,3,1,5,6] ; 
subplot(2,5,1) ; barwitherr(squeeze(std(mean(betas(:,:,stims),1),0,2))./sqrt(19),squeeze(mean(mean(betas(:,:,stims),1),2)))
subplot(2,5,2) ; barwitherr(squeeze(std(mean(hgs(:,:,stims),1),0,2))./sqrt(19),squeeze(mean(mean(hgs(:,:,stims),1),2)))

hz = 15:25 ; errorbarxy(squeeze(mean(mean(mtcomps(:,stims,hz),1),3)),squeeze(mean(mean(betas(:,:,stims),1),2)),squeeze(std(mean(mtcomps(:,stims,hz),3),0,1))./sqrt(16),squeeze(std(mean(betas(:,:,stims),1),0,2))./sqrt(16)) ;









%%% plot the %change masked by mean correlatoins
stims = [2,3,1,5,6] ; corrmask = squeeze(mean(corrs,1)) ;itvl = 60 ; icount = 1 ; x = 25:160 ; y = 20:200 ; itincr = 3 ;
for i=1:max(size(stims)) ; subplot(1,5,i) ; 
    mcorrs = squeeze(mean(corrs,1)) ; maskcorrs = mcorrs >.2 ;
    mperc = squeeze(mean(imgs(:,stims(i),:,:,:),1)) ; mperc = mperc.*(mcorrs>.15) ;
    img = fliplr(rot90(mat2gray(squeeze(mean(a(x,y,itvl:itvl+itincr),3)).^3))) ; 
    im2 = (squeeze(mean(mperc(x,y,itvl:itvl+itincr),3))) ; im2(1,1) = 0.015 ; 
    I = fliplr(rot90(uint8((mat2gray(im2)*255)))) ;  
    rgb = ind2rgb(gray2ind(I,255),jet(255)) ;
    imshow(img, 'InitialMag', 'fit') ;
    hold on ; h = imshow(rgb) ; hold off ; 
    set(h, 'AlphaData', fliplr(abs(rot90(squeeze(mean(mcorrs(x,y,itvl:itvl+itincr),3)))))>.15) ;
end

