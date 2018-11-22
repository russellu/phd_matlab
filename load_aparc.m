clear all ; close all ; 
%subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
%    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 

subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'} ; 


for sb=1:length(subs)
cd(['c:/shared/badger_struct/',subs{sb},'/txts']) ; 
txts=dir('*txt') ;
for txt=1:length(txts) 
fid = fopen([txts(txt).name]) ; 
labs = fgetl(fid) ; labsplit = strsplit(labs,' ') ; 
nums = fgetl(fid) ; numsplit = strsplit(nums,' ') ; 
for i=2:length(numsplit) ; vals(sb,txt,i-1) = str2num(numsplit{i}) ; end
for i=2:length(labsplit) ; labnames{txt,i-1} = labsplit{i} ; end
fclose(fid) ; 
end
end
small_vals = vals ; 
cd c:/shared/savedata_struct; save('small_vals','small_vals') ; save('labnames','labnames') ; 

