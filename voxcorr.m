function rvol = voxcorr(bold,timeseries)
% rvol = voxcorr(bold_image,timeseries) ; 
% function to correlate a 3d+t dataset (bold) with a 1d time series
% the time series must be the same size as the 4th dimension of time
% dimension of the BOLD image

corrvol = zeros(1,1,1,size(bold,4)) ; corrvol(1,1,1,:) = timeseries ;  
corrvol = repmat(corrvol,[size(bold,1),size(bold,2),size(bold,3),1]) ; 
meancorrvol = corrvol - repmat(mean(corrvol,4),[1,1,1,size(corrvol,4)]) ; % mean subtract 
meanbold = double(bold) - repmat(mean(bold,4),[1,1,1,size(bold,4)]) ;
% normalized co-variance
rvol = (sum(meancorrvol.*meanbold,4)./(sqrt(sum(meancorrvol.*meancorrvol,4)).*sqrt(sum(meanbold.*meanbold,4)))) ; 

end