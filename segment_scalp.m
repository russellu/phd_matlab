cd c:/shared/t1_ute ; ls 

t1 = load_untouch_nii('MONG_01_RB.nii.gz') ; 
t1img = t1.img ;
flat1 = reshape(t1img,[size(t1img,1),size(t1img,2)*size(t1img,3)]) ; 
[kc,km] = kmeans(round(mat2gray(flat1)*255),10) ; 
kimg = reshape(km,size(t1img)) ; 


[cx,cy,cz] = centmass3(kimg==1) ; 







