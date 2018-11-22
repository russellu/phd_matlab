function [labs,coords] = load_elocs(locname)
% function [labs,coords] = load_elocs
% load the electrode coordinate file output by the java program

cd c:/shared/elocs ; 
fid = fopen(locname) ; 
icount = 1 ; 
while ~feof(fid)
    line = fgets(fid) ;
    split = strsplit(line,',') ;
    labs{icount} = split{1} ; 
    coords(icount,:) = [str2num(split{2}),str2num(split{3})] ; 
    icount = icount + 1 ; 
end
fclose(fid)  ;
end

