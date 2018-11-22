%cd c:/users/butr2901/Documents/ ; 
cd C:\Users\butr2901\Documents\jforex
ls ; clear all ; close all ; 
ss = dir('savestats_*') ; 
for s = 1:length(ss)  ;
   a = load(ss(s).name) ; 
   [sv,si] = sort(a(:,5),'descend') ; 
   allvals(s,:,:) = a(si(1:500),:) ; 
   allnames{s} = strrep(ss(s).name,'savestats_','') ; 
end
%bar(squeeze(mean(allvals(:,:,5),2)))
%for i=1:size(allvals,1) ;figure; 
%    [sv,si] = (sort(allvals(i,:,5),'descend')) ; 
%    plot(allvals(i,si,1).*allvals(i,si,2)) ; hline(0,'k') ; title(allnames{i}); 
%end

% MAKE DIFFUSION FIGURE

cd c:/users/butr2901/documents/fxoutput30 ;
%sumconsist = squeeze((allvals(:,:,1).^2.*sqrt(allvals(:,:,2))).*allvals(:,:,5).^2) ; 
for i=1:size(allvals,1) ; 
    [sv,si] = sort(allvals(i,:,5),'descend') ;
    if sum(sv)>0
    dlmwrite([allnames{i},'.txt'],squeeze(allvals(i,si,6:end))) ; 
    figure,plot(squeeze(allvals(i,si,5))) ; 
    title(allnames{i}) ; 
    end
end