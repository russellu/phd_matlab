box = zeros(10,20) ; 
for i=1:size(box,2)
    if mod(i,2) == 0 
        box(:,i) = (i-1)*size(box,1)+1 : i*size(box,1) ; 
    else
        box(:,i) = fliplr((i-1)*size(box,1)+1 : i*size(box,1)) ; 
    end
end

 