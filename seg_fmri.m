cd c:/shared/allfmris/sub_jeremie ; 
f1 = load_untouch_nii('f1.nii.gz') ; 

[kc,km] = kmeans(f1.img(:),8,'MaxIter',500) ; 
kimg = reshape(kc,size(f1.img)) ; 
