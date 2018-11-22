clear all ; close all ; 
cd C:\shared\freesurfer_segs\sub_russell\SUMA ; ls 

names = {'bot_right','bot_left','top_right','top_left','left_hemi','right_hemi','top_hemi','bot_hemi','full','fov','periph'}; 
elecbrain = load_untouch_nii('elecbrain_T1.nii.gz'); 

for title=1:length(names)
   avg = load_untouch_nii(['avg_',names{title},'.nii.gz']); 
   allavgs(:,:,:,title) = avg.img; 
   
   [cx(title),cy(title),cz(title)] = centmass3(allavgs(:,:,:,title)>.2); 

end

for i=1:65
   [ex(i),ey(i),ez(i)] = centmass3(elecbrain.img==i);  
end

clear dists; 
for i=1:11
dists(i,:) = sqrt((cx(i)-ex).^2 + (cy(i)-ey).^2 + (cz(i)-ez).^2); 
end

cd e:/orientation_retinotopy/fatlas ; save('dists','dists'); 