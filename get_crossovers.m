function crossovers = get_crossovers(closedata,mvgdata)
% function get_crossovers(closedata,mvgdata)

crossovers = zeros(1,length(closedata)); 

for i=2:length(closedata)
   if closedata(i) > mvgdata(i) && closedata(i-1) < mvgdata(i-1)
       crossovers(i) = 1;
   elseif closedata(i) < mvgdata(i) && closedata(i-1) > mvgdata(i-1)
       crossovers(i) = -1; 
   end
end

end