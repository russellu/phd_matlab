function [Pi p ol] = Shepherd(a,b,j)
%
% [Pi p ol] = Shepherd(a,b,[j=10000])
%
% Calculates the robust Shepherd's Pi correlation:
%   It first bootstraps the Mahalanobis distances, removes all observations 
%   with m >= 6 and finally calculates Spearman's rho of the remaining data.
%
% The inputs a and b are column vectors containing the two variables.
% 
% The optional j sets the number of bootstraps (default=10000).  
%
% Pi is Spearman's Rho after outlier removal. 
% 
% p is multiplied by 2 to achieve a nominal false alarm rate. 
% The third output ol is a boolean vector indicating outliers.
%
% Requires the Statistics Toolbox.
% 

% Default number of bootstraps
if nargin < 3
    j = 10000;
end

% Bootstrapped Mahalanobis distances
if j > 0
    m = bsmahal([a b], [a b], j);
else
    m = mahal([a b], [a b]);
end

% Determine potential outliers
ol = (m >= 6);

% Calculate Pi
if sum(~ol) > 2
    [Pi p] = corr(a(~ol), b(~ol), 'type', 'spearman');
else
    Pi = NaN;
    p = NaN;
end

if j > 0
    % Adjust p
    p = p*2;
    if p > 1
        p = 1;
    end
end

