cd c:/shared/mnirefs ; ls 
bflows = dir('bloodflow*') ; 

for b=1:length(bflows)  ;
    nii = load_nii(bflows(b).name) ;  
    bimgs(b,:,:,:) = nii.img ; 
    
    
end