cd c:/shared/pcasl ; 

pcat = load_untouch_nii('mc_paslcat.nii.gz') ; 
pimg = pcat.img ; 
p1 = pimg(:,:,:,1:2:200) - pimg(:,:,:,2:2:200) ; 
p2 = pimg(:,:,:,201:2:400) - pimg(:,:,:,202:2:400) ; 
p3 = pimg(:,:,:,401:2:600) - pimg(:,:,:,402:2:600) ; 
allp(:,:,:,1:100) = p1 ; allp(:,:,:,101:200) = p2 ; allp(:,:,:,201:300) = p3 ; 
save_nii(make_nii(allp),'allp.nii.gz') 

meanp = (p1+p2+p3) /3 ; 
save_nii(make_nii(meanp),'meanp.nii.gz') ; 


pcat = load_untouch_nii('mc_pcasl.nii.gz') ; 
pimg = pcat.img ; 

p1 = pimg(:,:,:,1:2:128) - pimg(:,:,:,2:2:128) ; 
p2 = pimg(:,:,:,129:2:256) - pimg(:,:,:,130:2:256) ; 
p3 = pimg(:,:,:,257:2:384) - pimg(:,:,:,258:2:384) ; 
meanp = (p1 + p2 + p3) / 3 ;
save_nii(make_nii(meanp),'mean_pcasl.nii.gz') ; 








