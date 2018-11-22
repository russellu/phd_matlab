function emvg = getemvg(data,period)

npoints = period; incr = 1:5/npoints:5; 
kernel = exp(incr); kernel = kernel./sum(kernel); 
emvg = zeros(size(data)); 
for i=length(kernel)+1:length(data)
    emvg(i) = sum((data((i-length(kernel)+1):i)).*kernel');  
end


end