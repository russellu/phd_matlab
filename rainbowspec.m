%%% weight the frequencies according to color

T = pi ; 
n=50 ; 
fcount = 1 ; 
clear freqs ; 
for n=.1:.1:10
    freqs(fcount,:) = sin(1:n:1000*n) ; 
    
    fcount = fcount + 1 ; 
    
end


clear xs ;
xcount = 1 ; 
for x=1:.1:100
    xs(xcount,:) = sin(1:x:10*x) ;
    xcount = xcount + 1 ; 
end

clear xs
for i=1:10 ;
   xs(i,:) = rand(1,1000)*i ;   
end













