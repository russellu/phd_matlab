clear all ; close all ; 
cd c:/shared/Raw/subjects ; ls 
subs = dir('ica_S*set') ; 
for ss=1:length(subs) ; 
    EEG = pop_loadset(subs(ss).name) ; 
    figure ;
    for i=1:size(EEG.icawinv,2) ; subplot(5,13,i) ; topoplot(double(EEG.icawinv(:,i)),EEG.chanlocs,'electrodes','off') ; end ; title(subs(ss).name) ; 
end