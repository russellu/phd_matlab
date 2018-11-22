function ci = corrci(x, y, j, sr)
%
% ci = corrci(x, y, n, [j=10000, sr=false])
%
% Returns the nominal 95% confidence interval as well as the 
%   bootstrapped interval for Pearson's correlation between x and y.
%   
% j defines the number of bootstrap resamples (default = 10,000).
%
% If sr is true (default = false), the Spearman's rho correlation
%   coefficient is calculated instead of Pearson's r.
%
% To output is a 2x2 matrix: 
%   the first row contains the nominal confidence interval
%   the second row contains the bootstrapped interval
%
% Requires the Statistics Toolbox.
%

if nargin < 3
    j = 10000;
    sr = false;
elseif nargin < 4
    sr = false;
end

% Correlation coefficient
if sr
    r = corr(x,y,'type','spearman');
else
    r = corr(x,y);
end
n = size(x,1);

% Nominal interval
nci = tanh(atanh(r) + norminv(0.025)/sqrt(n-3)*[1 -1]);

% Bootstrapped interval
if sr
    Rs = bootstrp(j, @(x,y) corr(x,y,'type','spearman'), x, y);
else
    Rs = bootstrp(j, @(x,y) corr(x,y), x, y);
end
bci = prctile(Rs, [2.5 97.5]);

% Output both intervals
ci = [nci; bci];
