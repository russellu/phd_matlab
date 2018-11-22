function [xcent,ycent,zcent] = cmass3(img)
% get center of mass of 2d image
xs = 1:size(img,1) ; 
ys = 1:size(img,2) ; 
zs = 1:size(img,3) ; 

xvec = squeeze(sum(sum(img,2),3)) ; 
yvec = squeeze(sum(sum(img,1),3)) ; 
zvec = squeeze(sum(sum(img,1),2) ) ; %
%size(zvec)

xcent = round(sum(xvec.*xs')./sum(xvec)) ; 
ycent = round(sum(yvec.*ys)./sum(yvec)) ; 
zcent = round(sum(zvec.*zs')./sum(zvec)) ; 



end