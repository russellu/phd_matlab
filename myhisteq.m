
img = squeeze(t1im(:,:,65)) ; 
[x,zs] = myhist(img) ; 
csum = cumsum(x) ; 
nsum = round(((csum-min(csum))/(size(img,1)*size(img,2)-min(csum)))*255) ;
newhist = diff(nsum) ; 
newhist(size(newhist,2)+1) = 0 ; 


