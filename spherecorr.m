clear all ; close all ; 
sujets = {'sujet1','sujet2','sujet3','sujet4','sujet5','sujet6','sujet8','sujet9','sujet10','sujet11','sujet12'} ; 
for s=1:length(sujets) ; 
cd(['c:/shared/claudie/',sujets{s}]) ; ls 

postnii = load_untouch_nii('bp_post.nii.gz') ; 
prenii = load_untouch_nii('bp_pre.nii.gz') ; 
graynii = load_untouch_nii('3mm_cortex.nii.gz') ; 
motornii = load_untouch_nii('3mm_motor.nii.gz') ; 

cd spheres ; 
allspheres = dir('3mm_*') ; 
clear sphimgs
for sph=1:length(allspheres)
   sphnii = load_untouch_nii(['3mm_sphere_',num2str(sph),'.nii.gz']) ; 
   sphimgs(sph,:,:,:) = sphnii.img ;  
end

for sphind = 1:20%size(sphimgs,1)
    spherei = squeeze(sphimgs(sphind,:,:,:)) ; 
    [sphx,sphy,sphz] = ind2sub(size(spherei),find((spherei)>0)) ; %.*double(graynii.img)
    [handx,handy,handz] = ind2sub(size(spherei),find((double(motornii.img).*double(graynii.img))>0)) ; %
    clear pre_hand post_hand pre_sphere post_sphere
    for i=1:length(handx)
       pre_hand(i,:) = squeeze(prenii.img(handx(i),handy(i),handz(i),:)) ;  
       post_hand(i,:) = squeeze(postnii.img(handx(i),handy(i),handz(i),:)) ;  
    end
    for i=1:length(sphx)
       pre_sphere(i,:) = squeeze(prenii.img(sphx(i),sphy(i),sphz(i),:)) ;  
       post_sphere(i,:) = squeeze(postnii.img(sphx(i),sphy(i),sphz(i),:)) ;  
    end
    cpre = corr(pre_hand',pre_sphere') ; 
    cpost = corr(post_hand',post_sphere') ; 
    meanpre(sphind) = mean(mean(cpre)) ;
    meanpost(sphind) = mean(mean(cpost)) ; 
end

allpre(s,:) = meanpre ; 
allpost(s,:) = meanpost ; 


end

both = [mean(allpre(:,10:end),2),mean(allpost(:,10:end),2)] ; 


