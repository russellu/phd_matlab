function r = ndcorr(a,b)
% r = ndcorr(a,b) ; % where a,b are n-d arrays (n>=2)
% correlation coefficient (pearson) between two n-d arrays.
a = a(:) - mean(a(:)) ; b = b(:) - mean(b(:)) ; % mean subtract
r = sum(a.*b)./(sqrt(sum(a.*a)).*sqrt(sum(b.*b))) ; % normalized covariance
end
