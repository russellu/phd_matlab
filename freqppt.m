cd c:/users/acer/desktop ; ls
a = imread('fourier.png') ; 
a = rgb2gray(a) ; 
[x,y] = meshgrid(-size(a,2)/2:size(a,2)/2-1,-size(a,1)/2:size(a,1)/2-1) ; 
sig = 14 ; 
gaussf = exp(-(x.^2 + y.^2)/sig.^2) ;
igaussf = abs(x.*y) ; 
fa = fftshift(fft2(a)) ; 
fa = fa.*igaussf ; 
invfa = abs(ifft2(fa)) ; 
imagesc(invfa)