clear all ; close all; 
cd e:/switest/sorted/tr5s/mel ; ls 

mel = load_untouch_nii('melodic_IC.nii.gz'); 
mix = load('melodic_mix'); 
meanimg = load_untouch_nii('mean.nii.gz'); 

subplot(2,2,1); 
imagesc(rot90(squeeze(mean(meanimg.img(60:68,:,:),1)))); title('axial anat');
subplot(2,2,2); 
imagesc(rot90(squeeze(max(mel.img(55:73,:,:,10),[],1)))); colormap jet; title('axial artery component'); 

subplot(2,2,3); 
imagesc(rot90(squeeze(mean(meanimg.img(:,:,15:35),3)))); title('sagittal anat');
subplot(2,2,4); 
imagesc(rot90(squeeze(max(mel.img(:,:,15:35,10),[],3)))); colormap jet; title('sagittal artery component'); 

suptitle('melodic ICA results for 5s TR SWI'); 