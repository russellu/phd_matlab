clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','valerie','tegan'} ; 
for subj = 1:length(subs)
cd c:/shared/badger_eeg ; ls 
cd(subs{subj}) ; 
files = dir('*vhdr') ; bads = [] ; 
for i=1:length(files)
     if ~isempty(findstr('Pulse',files(i).name)) 
         bads(length(bads)+1) = i ; 
     end
end
files(bads) = []  ;
clear impvals ;
for f=1:length(files) ; 
fid = fopen(files(f).name) ;
a = textscan(fid,'%s') ; 
fstring = a{1} ; 
electrodes = cell(0) ; 
elecs = elecstring() ; 

for i=1:length(elecs)
    indi = find(strcmpi(elecs{i},fstring)) ; 
    if strcmpi(fstring(indi+1),'Out')
        disp('out of range') ; 
        elecvals(subj,f,i) = 100 ; 
    else 
        elecvals(subj,f,i) = str2double(fstring(indi+1)) ; 
    end
end

end

impedences = squeeze(elecvals(subj,:,:)) ; 
save('impedences','impedences') ; 

end


for i=1:9 ; subplot(3,4,i) ; imagesc(squeeze(elecvals(i,:,:)),[-100,100]) ; title(['mean impededance = ',num2str(mean(mean(elecvals(i,:,:)))),', ',subs{i}]) ; end
figure,barwitherr(squeeze(std(mean(elecvals,2),0,3))./sqrt(66),squeeze(mean(mean(elecvals,2),3))) ; 
% show topoplot of impedences




