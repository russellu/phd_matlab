function h = plotbrainRGB(anat,mask,heatmap,rotangle) 
% function h = plotoverlayIntensity2D(anat,mask,heatmap,rotangle) 

img = imrotate(uint8(mat2gray(anat)*255),rotangle) ; img(isnan(img)) = 0 ; img(isinf(img)) = 0 ; 
for i=1:3
    im2(:,:,i) = imrotate(uint8(mat2gray(heatmap(:,:,i))*255),rotangle) ; 
end
im2(isnan(im2)) = 0 ; im2(isinf(im2)) = 0 ; 
rgb = im2 ; imshow(rgb) ; 
imshow(img, 'InitialMag', 'fit') ;
hold on ; h = imshow(rgb) ; hold off ; 
set(h, 'AlphaData', imrotate((squeeze(mask)),rotangle)) ;

end