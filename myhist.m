function [h,zs] = myhist(img)

imin = min(min(img)) ; 
imax = max(max(img)) ;
h = zeros(1,imax-imin) ; 
zcount = 1 ;
zs = zeros(size(img,1),size(img,2),imax-imin) ; 
for i=imin:imax
    z = img==i ; 
    zs(:,:,zcount) = z ; 
    zcount = zcount + 1 ; 
    h(i-imin+1) = sum(sum(z)) ;   
end

end