function [retarr,datas] = get_dukascopy(path,id)
% function retarr = proc_dukascopy(path,id)
% cd to a directory containing dukascopy files get all files matching
% keyword 'id' and load them into an array

% move to directory and get all files 
cd(path) ;
datas = dir(id) ;

% loop over the files, put them into array
for i=1:size(datas)
    fid = fopen(datas(i).name) ;
    retarr{i} = textscan(fid,'%s%s%n%n%n%n%n') ;
    fclose(fid) ;
end

end