cd c:/users/acer/Documents/fx ; ls ; clear all ; close all ; 
ss = dir('savestats_*') ; 
for s = 1:length(ss)  ;
   a = load(ss(s).name) ;  
   allvals(s,:,:) = a ; 
   allnames{s} = strrep(ss(s).name,'savestats_','') ; 
end
% get the spreads
fid = fopen('c:/users/Acer/Documents/spreadTxt.txt') ; string = fgets(fid) ; spreadnames = {} ; spreads = [] ; 
while string ~= -1 ; 
    spliti = strsplit(string,',') ; spreadnames{length(spreadnames)+1} = spliti{1} ; 
    spreads(length(spreads)+1) = str2num(spliti{2}) ; string = fgets(fid) ; 
end
fclose(fid) ; 
% get the scalefactors
fid = fopen('c:/users/Acer/Documents/weightTxt.txt') ; string = fgets(fid) ; weightnames = {} ; weights = [] ; 
while string ~= -1 ; 
    spliti = strsplit(string,',') ; weightnames{length(weightnames)+1} = spliti{1} ; 
    weights(length(weights)+1) = str2num(spliti{2}) ; string = fgets(fid) ; 
end
fclose(fid) ; 

% cost function: (meanprofit-cost) > 0 .* (sumprofit.*consistency)










