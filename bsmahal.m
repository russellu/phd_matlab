function m = bsmahal(a,b,j)
%
% m = bsmahal(a,b,[j=10000])
%
% Bootstraps Mahalanobis distances of all pairs in a from the centre 
%   of the data cloud in b. Both a and b are 2-column matrices. 
%   j defines number of bootstrap samples to calculate (default = 10,000).
%
% Returns in m the Mahalanobis distance for each row in a, 
%   averaged across all the bootstrap resamples.
%
% Requires the Statistics Toolbox.
%
% (You could also use the MATLAB bootstrp function instead as this is 
%  much faster but this is buggy for large matrices so it is no good 
%  for calculating the contour plot in ScatterOutliers.m).
%

% Number of bootstraps
if nargin < 3
    j = 10000;
end

n = size(b,1);  % Number of data points
Ms = [];

% Bootstrap the distances
for i = 1:j
    x = ceil(rand(n,1)*n);
    s1 = b(x,1);
    s2 = b(x,2);
    m = mahal(a,[s1 s2]);
    Ms = [Ms m];
end

% Average across all bootstraps
m = mean(Ms,2);

