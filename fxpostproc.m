cd c:/users/butr2901/Documents ; ls 
ss = dir('savestats_*') ; 
for s = 1:length(ss)  ;
   a = load(ss(s).name) ;  
   [sv,si] = sort(a(:,5),'descend') ; 
   allmeans(s,:) = a(si,1) ; 
   allnts(s,:) = a(si,2) ; 
   allcosts(s,:) = a(si,5) ; 
   %figure,plot(a(si,1)) ; title(ss(s).name) ; 
   allnames{s} = strrep(ss(s).name,'savestats_','') ; 
end
spreadweights = zeros(1,length(allnames)) ; spreadweights(zscore(squeeze(mean(allmeans,2)))>0) = 0.01 ; spreadweights(zscore(squeeze(mean(allmeans,2)))<0) = 0.0001 ;
spreadweights(length(spreadweights)) = 0.001 ; 
fid = fopen('c:/users/butr2901/Documents/spreadFile.txt') ; 
fline1 = fgets(fid) ; fline2 = fgets(fid) ; fclose(fid) ; 
currs = strsplit(fline1,',') ; spreads = strsplit(fline2,',') ; 
currs = currs(1:length(currs)-1) ; spreads = (cellfun(@str2num,spreads(1:length(spreads)-1))) ; 


for i=1:length(allnames) ; ispreads(i) = spreads(find(strcmp(allnames{i},currs))) ; end

imagesc(squeeze(mean(allmeans,2))) ; set(gca,'YTick',1:44,'YTickLabel',allnames) ; 
figure,imagesc(zscore(squeeze(mean(allmeans,2)))>0) ; set(gca,'YTick',1:44,'YTickLabel',allnames) ; 
figure,imagesc(spreadweights), set(gca,'YTick',1:44,'YTickLabel',allnames) ; 
figure,imagesc(ispreads'.*spreadweights') ;  set(gca,'YTick',1:44,'YTickLabel',allnames) ; 

wspreads = ispreads.*spreadweights +0.7*spreadweights ; 
cmeans = allmeans - repmat(wspreads,[12500,1])' ; 
tps = cmeans.*allnts ; 
for i=1:44; 
    [sv,si] = (sort(tps(i,:),'descend')) ;
    figure,
    subplot(1,2,1) ; plot(tps(i,si)) ; hline(0,'k') ; 
    subplot(1,2,2) ; plot(allcosts(i,si)) ; 
    title(allnames{i}) ; 
end


