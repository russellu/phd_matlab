cd c:/shared/mnirefs2 ; ls 
corrsubs = dir('corr_*') ; 
for cs = 1:length(corrsubs)  ;
   nii = load_nii(corrsubs(cs).name) ; 
   cimgs(cs,:,:,:) = nii.img ; 
   
    
    
end

anii = load_nii('c:/shared/ATLASES/MNI152_T1_1mm.nii.gz') ;  
anat = anii.img ; 
meancorrs = squeeze(mean(cimgs,1)) ; 
plotoverlayIntensity(anat,mat2gray(meancorrs),meancorrs,70:75) ;