function [cx,cy] = centmass(img)
% get center of mass of 2d image
xs = (1:size(img,1))' ; size(xs)
ys = (1:size(img,2)) ; size(ys)


cx = floor(sum(sum(img,2).*xs)./sum(sum(img,2))) ;
cy = floor(sum(sum(img,1).*ys)./sum(sum(img,1))) ; 





end