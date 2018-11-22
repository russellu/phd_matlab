rands = rand(10,1) ; 
randinds = ceil(rand(10,1)*100) ; 
clear randmat
for i=1:size(rands,1)
    for j=1:size(rands,1)
        if j<i 
            randmat(i,j) = rands(i) - rands(j) ; 
        end
    end
end

% the path to the collapsed tree?
% give each group its own name, and save the clusters at each step. 














