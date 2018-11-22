cd c:/users/butr2901/pictures/; 
icon = imread('fxlogo.png'); 
newicon = zeros(512,512,3) ; 
for i=1:3 ; newicon(:,:,i) = imresize(icon(:,:,i),[512,512],'Method','nearest') ; end
newicon = uint8(mat2gray(newicon)*255) ;icon = newicon; 
%binicon = (icon>mean(icon(:))*.75); 
% imwrite(uint8(mat2gray(binicon)*255),'binicon.png'); 
alpha = double(ones(size(icon,1),size(icon,2))); 
meanicon = mean(icon,3); 
alpha(meanicon==255) = 0; 
alpha(meanicon==0) = 0.25; 
imwrite(icon,'trans_forex_logo2.png','png','Alpha',alpha); 