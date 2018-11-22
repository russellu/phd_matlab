function padded = pad3d(input,padamt) 
% function padded = pad3d(input,padamt) ; 
% pad a 3d input array with zeros
padded = zeros(size(input)+padamt*2) ; 
padded(padamt:size(input,1)+padamt-1,padamt:size(input,2)+padamt-1,padamt:size(input,3)+padamt-1) = input ; 
end