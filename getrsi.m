% function rsi = getrsi(mopen,mclose,period)
function rsi = getrsi(mopen,mclose,period)

rsi = zeros(size(mopen)); 
for i=period+1:length(mopen)
    meangain = 0;
    meanloss = 0;  
    for j=i-period:i
        pricediff = mclose(j)-mopen(j);       
        if pricediff>0
            meangain = meangain + pricediff;
        elseif pricediff < 0
            meanloss = meanloss - pricediff; 
        end
    end
    meangain = meangain / period;
    meanloss = meanloss / period; 
    rsi(i) = (100 - 100/(1+(meangain/meanloss))); 
end


end



