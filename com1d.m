function cx = com1d(vec)
% get the center of mass of a 1d vector

if size(vec,1) == 1
    sz = 2 ; 
else sz = 1 ; 
end

if sz ==1 
cx = sum(vec'.*(1:size(vec,sz)))./sum(vec) ; 
else
    cx = sum(vec.*(1:size(vec,sz)))./sum(vec) ; 
end

end