function padded = zeropad(arr,amt)
% zero pad an array
% function padded = zeropad(arr,amt)


if ndims(arr) == 2
    paddy = zeros(size(arr,1)+amt*2,size(arr,2)+amt*2) ; 
    paddy(amt+1:size(paddy,1)-amt,amt+1:size(paddy,2)-amt) = arr ; 
elseif ndims(arr) == 3 
   for i=1:size(arr,3) % assume third dimension is color channel 
        paddy(:,:,i) = zeros(size(arr,1)+amt*2,size(arr,2)+amt*2) ; 
        paddy(amt+1:size(paddy,1)-amt,amt+1:size(paddy,2)-amt,i) = arr(:,:,i) ; 
   end
end

padded = paddy ;


end