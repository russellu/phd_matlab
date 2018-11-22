clear all ; close all ; 
cd c:/users/butr2901/documents/savestats ; ls 
ss = dir('savestats*') ; 
for s=1:length(ss)
    allss(s,:,:) = load(ss(s).name) ; 
    names{s} = strrep(ss(s).name,'savestats_','') ; 
    
end

imagesc(squeeze(allss(:,:,4)),[.5,1]) ; set(gca,'YTick',1:length(names),'YTickLabel',names) ; 