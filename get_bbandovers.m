function bbandover = get_bbandover(closedata,bbands)
%function bbandcrossover = get_bbandcrossovers(closedata,bbands)

bbandover = zeros(1,length(closedata)); 

for i=2:length(closedata)    
    if closedata(i) > bbands(1,i) 
        bbandover(i) = 1; 
    elseif closedata(i) < bbands(2,i)
        bbandover(i) = -1; 
    end   
end

end