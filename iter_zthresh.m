function [newdata,bads] = iter_zthresh(input,zthresh) 
% perform interative zscore thresholding
% receives a vector (size == (1,n)) and iteratively removes outliers. 
% zthresh is the absolute value of the allowed input

bads = zeros(size(input)) ; 

while max(zscore(input)) > zthresh %|| min(zscore(input)) < -zthresh ; 
    meaninput = mean(input) ; 
    zinput = zscore(input) ; 
    badies =  zinput > zthresh ;%| zinput < -zthresh  ; 
    bads = bads + badies ; 
    goodies = ~badies ; % good is not bad
    interpinput = input ;
    interpinput(badies) = meaninput ;
    input = interpinput ; 
end

newdata = input ; 



end