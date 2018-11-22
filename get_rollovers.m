cd C:\shared
[x,y,full] = xlsread('swaps.xlsx'); 
shorts = cell2mat({full{:,3}}); 
longs = cell2mat({full{:,2}}); 

