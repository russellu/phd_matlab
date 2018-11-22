[xg,yg] = meshgrid(-200:200) ; 
% make the moving bar 


for i=-200:200
    width = sqrt(abs(i))+5 ; 
    bar = xg>i & xg <i+width ;
    imagesc(bar) ;
    getframe ; 
   
end









