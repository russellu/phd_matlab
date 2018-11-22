cd C:\shared\epireg\rsns
corrimgs = load_untouch_nii('corrimgs.nii.gz') ;

for i=1:100 ; 
    figure,
    imgi = squeeze(corrimgs.img(:,:,:,i)); 
    subplot(2,2,1);
    imagesc(rot90(squeeze(max(imgi,[],1))));    
    subplot(2,2,2);
    imagesc(rot90(squeeze(max(imgi,[],2)))); 
    subplot(2,2,3); 
    imagesc(imrotate(squeeze(max(imgi,[],3)),270)); 
end
names = {'higher_visual','left_executive','anterior_visual','auditory','anterior_dmn_lobe','right_executive','motor','dmn','primary_visual','precuneus',...
    'language','somatosensory'};
indices = [44,41,38,33,31,30,28,15,13,96,85,76]; 
meanepi = load_untouch_nii('c:/shared/epireg/mean.nii.gz'); 
mkdir networks; cd networks; 
for n=1:length(names)
    imgi = corrimgs.img(:,:,:,indices(n)); 
    meanepi.img = imgi; 
    save_untouch_nii(meanepi,names{n}); 
end





