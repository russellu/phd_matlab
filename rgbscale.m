function scaled = rgbscale(img,sig)
% function scaled = rgbscale(img) ;
% splits a 2d+T img into an rgb img, turns the t whatever it might be into
% only 3 channels by linear weighting of indices

g = -1:2/size(img,3):1 ; 
g= exp(-(g.^2)./sig) ;
g = imresize(g,[1,size(img,3)]) ; 
r = g((size(g,2)/2):size(g,2)) ; 
r = imresize(r,[1,size(img,3)]) ;
b = g(1:size(g,2)/2) ; 
b = imresize(b,[1,size(img,3)]) ; 

redIM = zeros(1,1,size(img,3)) ; redIM(1,1,:) = r ; 
greenIM = zeros(1,1,size(img,3)) ; greenIM(1,1,:) = g ; 
blueIM = zeros(1,1,size(img,3)) ; blueIM(1,1,:) = b ;

redIM = repmat(redIM,[size(img,1),size(img,2),1]) ; 
greenIM = repmat(greenIM,[size(img,1),size(img,2),1]) ; 
blueIM = repmat(blueIM,[size(img,1),size(img,2),1]) ; 

scaled(:,:,1) = sum(mat2gray(redIM.*img),3) ; 
scaled(:,:,2) = sum(mat2gray(greenIM.*img),3) ; 
scaled(:,:,3) = sum(mat2gray(blueIM.*img),3) ; 

end