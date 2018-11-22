clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','valerie'} ;

cd(['c:/shared/badger_eeg/',subs{7}]) ; ls 
preprocs = dir('preproc*set') ; 
for i=1:length(preprocs)
EEG = pop_loadset(preprocs(i).name) ; 
if i==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged,1) ; end 
end

[s,f] = spectopo(merged.icaact,0,merged.srate,'plot','off') ; 
figure,
for i=1:64 ; 
    subplot(5,13,i) ;
    plot(squeeze(s(i,1:200)),'LineWidth',2) ; xlim([1,200]) ; title(i) ; ylim([-30,20]) ; 
    maxind = find(s(i,f>5 & f<20)==max(s(i,f>5 & f<20)),1) ; 
    vline(length(find(f<5))+maxind,'k') ; peaks(i) = length(find(f<5))+maxind ; 
end
figure,for i=1:64 ; subplot(5,13,i) ; topoplot(merged.icawinv(:,i),merged.chanlocs) ; title(f(peaks(i))) ; end

