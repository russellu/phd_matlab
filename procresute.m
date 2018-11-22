cd c:/shared/fute ; ls 
masked = load_untouch_nii('notnotmask.nii.gz') ; 
notmask = ~masked.img ; 
bw = bwconncomp(notmask) ; 
bwlist = bw.PixelIdxList ; 
zmask = zeros(size(notmask)) ; 
zmask(bwlist{1}) = 1 ; 
masked.img = zmask ; save_untouch_nii(masked,'notnotmask2.nii.gz') ; 

