cd c:/shared/utef/jeremie; ls 

locs = load_untouch_nii('color_locs.nii.gz') ; 
raw = load_untouch_nii('dof12_jeremie.nii.gz') ; 

for i=1:65
    [cx,cy,cz] = centmass3(locs.img==i) ; 
    locimg(i,:,:,:) = squeeze(raw.img(cx-10:cx+10,cy-10:cy+10,cz-10:cz+10)) ; 
end


for i=1:65 ; figure ; for j=1:21 ; subplot(4,6,j) ; imagesc(squeeze(locimg(i,j,:,:))) ; end ; end
