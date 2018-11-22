function mvg = getmvg(data,period)
%function mvg = getmvg(data,period)

for i=period+1:length(data)
    if i==period+1
        mvg(i) = mean(data(i-period+1:i)) ; 
    else
        mvg(i) = mvg(i-1) + data(i)/(period) - data(i-period)/period ;
    end  
end

end