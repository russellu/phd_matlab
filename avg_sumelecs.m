clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','tegan','valerie'} ;
for sb=1:length(subs)
    cd(['c:/shared/all_white_normals/a2_good/sub1_',subs{sb}]) ;    
    coords = load_untouch_nii('std_dilcoords.nii.gz') ; 
    length(unique(coords.img(:)))
    for i=1:65
       [cx(sb,i),cy(sb,i),cz(sb,i)] = centmass3(coords.img==i) ;  
    end
    
    
end
mcx = round(mean(cx,1)) ; mcy = round(mean(cy,1)) ; mcz = round(mean(cz,1)) ; 
sumsumreg = load_untouch_nii('c:/shared/ATLASES/sumsumreg.nii.gz') ; 
zs = zeros(size(sumsumreg.img)) ; 
for i=1:length(mcx) ; zs(mcx(i)-2:mcx(i)+2,mcy(i)-2:mcy(i)+2,mcz(i)-2:mcz(i)+2) = i ; end 
sumsumreg.img = zs ; cd c:/shared/ATLASES ; save_untouch_nii(sumsumreg,'sumsumreg_elecs.nii.gz') ; 

