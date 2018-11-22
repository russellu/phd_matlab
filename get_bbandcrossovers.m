function bbandcrossover = get_bbandcrossovers(closedata,bbands)
%function bbandcrossover = get_bbandcrossovers(closedata,bbands)

bbandcrossover = zeros(1,length(closedata)); 

for i=2:length(closedata)    
    if closedata(i) > bbands(1,i) && closedata(i-1) < bbands(1,i-1)
        bbandcrossover(i) = 1; 
    elseif closedata(i) < bbands(2,i) && closedata(i-1) > bbands(2,i-1)
        bbandcrossover(i) = -1; 
    end   
end

end