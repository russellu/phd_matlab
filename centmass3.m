% function [cx,cy,cz] = centmass3(img) 
% returns the 3d center of mass of the volume
function [cx,cy,cz] = centmass3(a)
x = 1:size(a,1) ; sumx = sum(sum(a,2),3).*x' ; cx = floor(sum(sumx./sum(sum(sum(a))))) ;
y = 1:size(a,2) ; sumy = sum(sum(a,1),3).*y ; cy = floor(sum(sumy./sum(sum(sum(a))))) ;
z = 1:size(a,3) ; sumz = squeeze(sum(sum(a,1),2)).*z' ; cz = floor(sum(sumz./sum(sum(sum(a))))) ;
end