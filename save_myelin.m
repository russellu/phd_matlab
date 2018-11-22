save_nii(make_nii(maps.gdn),'myelin_gdn.nii.gz') ; 
save_nii(make_nii(maps.ggm),'myelin_ggm.nii.gz') ; 
save_nii(make_nii(maps.gva),'myelin_gva.nii.gz') ; 
save_nii(make_nii(maps.alpha),'myelin_alpha.nii.gz') ; 
save_nii(make_nii(maps.FNR),'myelin_FNR.nii.gz') ; 


subplot(2,2,1) ; imagesc(maps2.sfr(:,:,12),[0,0.2]) ; title('sfr')  ; 
subplot(2,2,2) ; imagesc(maps2.sgm(:,:,12)) ; title('sgm')  ; 
subplot(2,2,3) ; imagesc(maps2.mfr(:,:,12)) ; title('mfr')  ; 
subplot(2,2,4) ; imagesc(maps2.mgm(:,:,12)) ; title('mgm')  ; 


for i=1:40  ;
    subplot(4,10,i) ; 
    imagesc(smoothmyelin(:,:,i),[0,0.2]) ;
end


myelin = maps2.sfr ; smoothmyelin = medfilt3(myelin) ; 