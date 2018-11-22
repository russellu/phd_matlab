cd C:\shared\raw\MONG_01_RB\ica_5
anat = load_untouch_nii('to_f1_f_mc_1.nii.gz') ; anatimg = anat.img ; 
[cx,cy,cz] = centmass3(anatimg) ; 
ica = load_untouch_nii('melodic_IC.nii.gz') ; icaimg = ica.img ; 
ts = load('melodic_mix') ; 
for i=1:size(ts,2)
   figure,
   subplot(2,2,1) ; imagesc(squeeze(anatimg(:,:,cz))) ; 
   plotoverlayIntensity2D(squeeze(anatimg(:,:,cz)),squeeze(mean(icaimg(:,:,cz-5:cz+5,i),3)),squeeze(mean(icaimg(:,:,cz-5:cz+5,i),3)),0) ; 
   subplot(2,2,2) ;    plotoverlayIntensity2D(squeeze(anatimg(cx,:,:)),squeeze(mean(icaimg(cx-5:cx+5,:,:,i),1)),squeeze(mean(icaimg(cx-5:cx+5,:,:,i),1)),90) ;  
   subplot(2,1,2) ; plot(squeeze(ts(:,i))) ; 
   
   suptitle(num2str(i)) ;  
end






