cd E:\switest

tof = load_untouch_nii('bet.nii.gz'); 

subplot(1,3,1) ;imagesc(max(tof.img,[],3)); colormap jet 
subplot(1,3,2) ;imagesc(rot90(squeeze(max(tof.img,[],2)))); 
subplot(1,3,3) ;imagesc(rot90(squeeze(max(tof.img,[],1)))); 

