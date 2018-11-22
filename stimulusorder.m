%%%% stimulus order: keep the stimuli apart from each other by atleast one
%%%% presentation (or two, if there are enough parameters)



params = 1:10 ; 
npresents = 4 ; 
randparams = [1,1,1,1] ; 
while (min(abs(diff(randparams))) < 1) | (min(abs(diff(randparams(1:2:end)))) < 1) | (min(abs(diff(randparams(2:2:end)))) < 1)
paramrep = repmat(params,[1,npresents]) ; 
randinds = randperm(length(paramrep)) ; 
randparams = paramrep(randinds) ; 
end








