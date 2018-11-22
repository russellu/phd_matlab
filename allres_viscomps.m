clear all ; close all 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','lisa','marc','marie','mathieu',...
    'maxime','mingham','patricia','po','russell','sunachakan','tah','thititip','vincent'} ; 

comps = {[18,13],[13,10],[10,8],[13,7],[23,12],[15,6],[18,8],[27,30],[16,20],[18,9],[10,9],[29,15],[6,8],[10,12],...
    [14,36],[13,8],[19,34],[8,16],[21,20],[9,13],[23,16],[17,12],[18,7],[11,25]} ;

for sub=24 ; 
cd(['c:/shared/all_white_normals/fmris/sub_',subs{sub},'/melodic/']) ; 
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
    subplottight(2,1,2) ;  plot(squeeze(mix(1:end,comp))) ; 
    suptitle(['comp',num2str(comp)]) ; 
end

end


