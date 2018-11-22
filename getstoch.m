function stoch = getstoch(data,period)
%function stoch = getstoch(data,period)

stochper = period; stoch = zeros(size(data)); 
for i=stochper+1:length(data)
    high = max(data(i-stochper:i)); 
    low = min(data(i-stochper:i)); 
    stoch(i) = 100*((data(i)-low)/(high-low));
end

end