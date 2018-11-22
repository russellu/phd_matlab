function wmvg = getwmvg(data,period)

weights=1:period; weights = weights./sum(weights); 
wmvg = zeros(size(data)); 
for i=period+1:length(data)
    datai = data(i-period+1:i); 
    wmvg(i) = sum(datai.*weights'); 
end

end
