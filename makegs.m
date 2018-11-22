cd c:/shared/lastute ; 
subs=dir('*') ; 
subs(1:2) = [] ; 
for s=4:length(subs)
   cd(['c:/shared/lastute/',subs(s).name]) ;  
   resute = load_untouch_nii('res_ute.nii.gz') ; 
   gs = resute.img - imfilter(resute.img,fspecial('gaussian',20,10)) ; 
   resute.img = gs ; 
   save_untouch_nii(resute,'gs_20.nii.gz') ; 
    
    
    
    
end

% avg
%{
cd c:/shared/regute2 ; ls
resute = load_untouch_nii('sumsumreg.nii.gz') ; 
gs = resute.img - imfilter(resute.img,fspecial('gaussian',20,10)) ; 
resute.img = gs ; 
save_untouch_nii(resute,'refute_20.nii.gz') ; 
%}






