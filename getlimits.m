function [xl, yl] = getlimits(xvec, yvec)
% function [xl, yl] = getlimits(xvec, yvec)

xl = [min(xvec),max(xvec)];
yl = [min(yvec),max(yvec)];

xdiff = (xl(2) - xl(1))/5; 
ydiff = (yl(2) - yl(1))/5; 

xl(1) = xl(1) - xdiff; 
xl(2) = xl(2) + xdiff;
yl(1) = yl(1) - ydiff;
yl(2) = yl(2) + ydiff; 

end