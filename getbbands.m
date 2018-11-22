function bbands = getbbands(data,period,stdfactor)
% function bbands = getbbands(data,period,stdfactor)

mvg = getmvg(data,period); 
bbands = zeros(2,length(data)); 
for i=period+1:length(data)
    stdi = std(data(i-period:i));
    bbands(1,i) = mvg(i)+stdi*stdfactor;
    bbands(2,i) = mvg(i)-stdi*stdfactor; 
end

end