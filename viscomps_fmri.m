cd(['c:/shared/allfmris/sub_gina/mel']) ;
all=dir('*') ; all(1:2) = [] ; 
meanNII = load_untouch_nii('mean.nii.gz') ; 
[cx,cy,cz] = centmass3(meanNII.img) ;   
comps = load_untouch_nii('melodic_IC.nii.gz') ;     
axial_anat = squeeze(mean(meanNII.img(:,:,cz-5:cz+5),3)) ; 
sagittal_anat = rot90(squeeze(mean(meanNII.img(cx-5:cx+5,:,:),1))) ; 
mix = load('melodic_mix') ; 
for comp = 1:size(comps.img,4)
    figure ; 
    subplottight(2,2,1) ;
    axial = (mean(squeeze(comps.img(:,:,:,comp)),3)) ; 
    plotoverlayIntensity2D(axial_anat,sqrt(mat2gray(axial)),axial,0) ; 
    subplottight(2,2,2) ; 
    sagittal = rot90(squeeze(mean(squeeze(comps.img(:,:,:,comp)),1))) ; 
    plotoverlayIntensity2D(sagittal_anat,sqrt(mat2gray(sagittal)),sagittal,0) ; 
    subplottight(2,1,2) ;  plot(squeeze(mix([1:245],comp))) ; 
    suptitle(['comp',num2str(comp)]) ; 
end
%{
goodcs = [45,30,29,] ;
badcs = zeros(1,size(comps.img,4)) ; badcs(goodcs) = 1 ; badcs = find(badcs==0) ; 
dlmwrite('badcs.txt',badcs) ; 
%}



