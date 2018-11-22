% function cci = getcci(low,high,close,period)

function cci = getcci(low,high,close,period)

    typrice = (low + high + close) / 3; 
    ccimvg = getmvg(typrice,period); 
    cci = zeros(size(low)); 
    
    for i=period+1:length(low)
        madi = mad(typrice(i-period:i)); 
        cci(i) = (1/0.015)*((typrice(i)-ccimvg(i))/madi); 
    end

end