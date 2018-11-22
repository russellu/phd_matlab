cd c:/shared/freesurfer_segs/sub_russell/SUMA ; ls 

lowcards = load_untouch_nii('low_cardinal.nii.gz'); 
highobs = load_untouch_nii('high_oblique.nii.gz') ;


newimg = zeros(size(highobs.img)); 
newimg(highobs.img>2.5) = 2;
newimg(lowcards.img>2.5) = 1;

lowcards.img = newimg ;
save_untouch_nii(lowcards,'color_oblique_card.nii.gz'); 