% function aroonosc = getaroon(data,aroonperiod)
function aroonosc = getaroon(data,aroonperiod)

aroonvals = zeros(2,length(data)); 
for i=aroonperiod+1:length(data)
    maxind = find(data(i-aroonperiod:i)==max(data(i-aroonperiod:i)),1); 
    minind = find(data(i-aroonperiod:i)==min(data(i-aroonperiod:i)),1); 
    aroonvals(1,i) = (maxind/aroonperiod)*100;
    aroonvals(2,i) = (minind/aroonperiod)*100; 
end
aroonosc = aroonvals(1,:) - aroonvals(2,:); 

end