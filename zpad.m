function padded = zpad(input,amt)
% function padded = zpad(input,amt)
% symmetrically pad an array with amt on all dimensions.

if ndims(input) == 3 % 3d case
    padded = zeros(size(input,1)+amt*2,size(input,2)+amt*2,size(input,3)+amt*2) ; 
    padded(amt:size(input,1)+amt-1,amt:size(input,2)+amt-1,amt:size(input,3)+amt-1) = input ; 
end

end