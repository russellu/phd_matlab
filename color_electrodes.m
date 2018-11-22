cd c:/shared/raw/MONG_01_RB/
gs = load_untouch_nii('gs.nii.gz') ; gsimg = gs.img ; 
nlayers = load_untouch_nii('native_layers.nii.gz') ; 
skinter = load_untouch_nii('skinter.nii.gz') ; 
skin = load_untouch_nii('skin.nii.gz') ; 
elabels = load('elabels') ; elabels = elabels.elabels ;
cskinter = load_untouch_nii('cskinter.nii.gz') ; cskimg = cskinter.img ; 
gsimg = gsimg.*imerode(nlayers.img > 0,strel(ones(3,3,3))) ;
wsize = 8 ; clear echunks ; indimg = zeros(size(gsimg)) ;
for i=1:63 ; 
   [cx,cy,cz] = centmass3(cskimg==i) ; 
   echunks(i,:,:,:) = gsimg(cx-wsize:cx+wsize,cy-wsize:cy+wsize,cz-wsize:cz+wsize) ;
   indimg(cx-wsize:cx+wsize,cy-wsize:cy+wsize,cz-wsize:cz+wsize) = i ; 
   chunkinds(i,:) = find(indimg==i) ; 
end

elecimg = zeros(size(gsimg)) ; 
for i=1:63
    vx = squeeze(echunks(i,:,:,:)) ; 
    vxz = medfilt3(vx) ; vxz = (vxz > mean(mean(mean(vxz)))*4).*(vxz~=0) ; 
    %allvx = vx(vx~=0) ; allvxinds = find(vx~=0) ; 
    %vxz(allvxinds(zscore(allvx)>3)) = 0 ;
    elecimg(chunkinds(i,:)) = vxz(:) ; 
end
gs.img = elecimg ; save_untouch_nii(gs,'elecimg.nii.gz') ; 




