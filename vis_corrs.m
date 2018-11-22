cd c:/shared/epi_atlas ; ls 
meanrest = load_untouch_nii('mm_mean_single2.nii.gz') ;

anat = load_untouch_nii('c:/shared/epi_atlas/ss.nii.gz') ; 

beta = squeeze(mean(meanrest.img(:,:,:,40:90),4)) ; 
beta(anat.img==0) = 0 ; 
beta(1,1,:) = -.2 ; beta(1,2,:) = .2 ; 
icount =1 ;
for i=50:4:70
    subplottight(1,6,icount) ; 
    plotoverlayIntensity2D(squeeze(anat.img(:,:,i)),mat2gray(abs(squeeze(beta(:,:,i)))),squeeze(beta(:,:,i)),270) ; icount = icount +1 ; 
end

imagesc([-.2,.2]) ; colorbar 