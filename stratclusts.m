%%%clustering by stratification based on variance within and between
%%%clusters

rvec = rand(10,1) ; 
nclusts = 3 ; 
% find all the combinations
% constraint: each element must be in a cluster, the same element cannot
% repeat twice. 