cd('c:/shared/raw/MONG_01_RB') ; 
layers = load_untouch_nii('native_layers.nii.gz')  ;
layers = layers.img ; 
rgb(:,:,:,1) = layers<=2 & layers>0 ; rgb(:,:,:,2) = layers>2 & layers <=4 ; rgb(:,:,:,3) = layers>4 & layers<=6 ; 
rgb = uint8(mat2gray(rgb)*255) ; 

gs = load_untouch_nii('gs.nii.gz') ; 
gsimg = gs.img-imfilter(gs.img,fspecial('gaussian',31,31)) ; 

rotangle = 90 ; 
anat = squeeze(gsimg(120,:,:)) ; 
rgbimg = squeeze(rgb(120,:,:,:)) ; for i=1:3 ; newrgb(:,:,i) = imrotate(squeeze(rgbimg(:,:,i)),rotangle) ; end
sumrgbimg = squeeze(sum(rgbimg,3)) ; 
mask = sumrgbimg>9999 ; 
img = imrotate(uint8(mat2gray(anat)*255),rotangle) ; img(isnan(img)) = 0 ; img(isinf(img)) = 0 ; 
%im2 = imrotate(uint8(mat2gray(heatmap)*255),rotangle) ; im2(isnan(im2)) = 0 ; im2(isinf(im2)) = 0 ;
%I = uint8(mat2gray(im2)*255) ;  
%rgb = ind2rgb(gray2ind(I,255),jet(255)) ; imshow(rgb) ;
imshow(newrgb) ; 
imshow(img, 'InitialMag', 'fit') ;
hold on ; h = imshow(newrgb) ; hold off ; 
set(h, 'AlphaData', imrotate((squeeze(mask)),rotangle)) ;


