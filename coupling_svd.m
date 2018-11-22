clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
for sub=9%:length(subs)
cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls
mov = dir('bcgica*rest*set') ; 
EEG = pop_loadset(mov(1).name) ; 
EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
rawts = mean(EEG.data([61,29],:),1) - mean(EEG.data([62,30],:),1) ;
chan = rawts ;

corrs = corr(abs(rawts)',abs(EEG.data)') ; 

x = cell2mat({EEG.chanlocs.X}) ; y = cell2mat({EEG.chanlocs.Y}) ; z = cell2mat({EEG.chanlocs.Z}) ; 
locs = [x;y;z] ;
blocs = zeros(3,64) ; blocs(:,1:31) = locs(:,1:31) ; blocs(:,33:end) = locs(:,32:end) ; 
diffs = zeros(64,64)  ;
for i=1:size(blocs,2)
    for j=1:size(blocs,2)
        diffs(i,j) = sqrt(sum((blocs(:,i)-blocs(:,j)).^2)) ; 
    end
end
for i=1:64 ; subplot(5,13,i) ; topoplot(diffs(i,:),EEG.chanlocs,'electrodes','numbers') ; title(i) ; end 



smth = (smooth(abs(rawts),500)) ;
plot(smth) ; hline(median(smth)*1.5) ; hold on ; plot(abs(rawts),'k') ; 
badinds = find(smth>median(smth)*1.5) ; 
bdatas = EEG.data(:,badinds) ; 

corrbdatas = double(corr(bdatas')) ; 
figure,for i=1:64 ; subplot(5,13,i) ; topoplot(corrbdatas(:,i),EEG.chanlocs,'maplimits',[-1,1]) ; end
thresh = 15 ; 
subinds = find(sum(corrbdatas.^2)>thresh) ; 
for i=1:length(subinds)
    
    
end




%[weights,sphere] = runica(bdatas) ; 
%winv = pinv(weights*sphere) ; for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),EEG.chanlocs) ; end


ica = pop_runica(EEG,'runica') ; 





end

freqs = 1:80 ;
allpows = zeros(size(ica.icaact,1),length(freqs),size(ica.icaact,2)) ; 
for f=1:length(freqs)
    allpows(:,f,:) = eegfiltfft(ica.icaact,250,freqs(f)-2,freqs(f)+2) ;   
end




