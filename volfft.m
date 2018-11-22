cd c:/shared
t1 = load_untouch_nii('Test_fMRI_2014_09_19_WIP_sT1W_3D_TFE_tra_SENSE_11_1.nii') ;
t1im = double(t1.img) ;
t1f = fftshift(fftn(t1im))  ;
[cx,cy,cz] = ndgrid(-size(t1im,1)/2:size(t1im,1)/2-1,-size(t1im,2)/2:size(t1im,2)/2-1,-size(t1im,3)/2:size(t1im,3)/2-1) ;
clear shells
shellsize = 10 ; 
shellcount = 1 ;
for i=10:5:75
    shells(shellcount,:,:,:) = sqrt(cx.^2 + cy.^2 + cz.^2) > i & sqrt(cx.^2 + cy.^2 + cz.^2) < i+shellsize ; 
    sig = i ; 
    gaussmask = exp(-(cx.^2+cy.^2+cz.^2)/sig) ; 
    t1fs(shellcount,:,:,:) = t1f.*gaussmask ; 
    shellcount = shellcount + 1 ;  
    disp(shellcount) ;
end
rawshells = round(abs(ifftn(squeeze(t1fs(2,:,:,:))))) ; 
clear shellthreshs
tcount = 1 ;
for i=1:5:116 ; 
    shellthreshs(tcount,:,:,:) = rawshells > i ;
    a = make_nii(double(squeeze(shellthreshs(tcount,:,:,:)))) ;
    save_nii(a,[num2str(i),'.nii.gz']) ; 
    tcount = tcount + 1 ; 
end
a = make_nii(t1im) ; 
save_nii(a,'0a.nii.gz') ; 

