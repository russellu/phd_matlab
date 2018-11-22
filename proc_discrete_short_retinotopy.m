clear all ; close all ; 
cd C:\shared\discret
discr = dir('disc*vhdr') ; 
for i=1:length(discr)
   EEG = pop_loadbv('.',discr(i).name) ; 
   EEG = pop_resample(EEG,256) ; 
   if i==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG) ; end
end
merged = pop_chanedit(merged,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
temp = merged.data(1:32,:) ; merged.data(1:32,:) = merged.data(33:64,:) ; merged.data(33:64,:) = temp ; 
mergefilt = merged ; mergefilt.data = eegfiltfft(merged.data,merged.srate,1,90) ; 
ica = pop_runica(mergefilt,'runica') ; ica.icaact = ica.icaweights*ica.icasphere*mergefilt.data ; 
ica.data = ica.icaact ; 
stims = {'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8','S  9','S 10','S 11','S 12','S 13','S 14','S 15','S 16','S 17','S 18','S 19','S 20','S 21','S 22','S 23','S 24','S 25','S 26','S 27','S 28'} ;
clear allersp ; 
for s=1:length(stims) ; disp(s) ; 
    allep = pop_epoch(ica,{stims{s}},[-.5,2]) ;
    for i=1:64 ;
            [allersp(s,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off') ; 
                
    end
end
for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(allersp(:,i,:,:),1)),[-3,3]) ; title(i) ; end


goodcs = [6,12,14,17,23,27,29,32] ; 
mtersp = squeeze(mean(mean(allersp(:,goodcs,:,times>0.25 & times<1),2),4)) ; 

bads = zeros(1,64) ; bads(goodcs) = 1 ; bads = find(bads==0) ; 
acts = ica.icaact ; acts(bads,:) = 0 ; 
invacts = pinv(ica.icaweights*ica.icasphere)*acts ; 
ica.data = invacts ; postelecs = [60,29,30,31,64,63,62,61,23,56,24,57,25,58,26,59,27,28,32] ; 
freqs = 1:10:100 ; alldats = zeros(28,64,length(freqs),333,90) ; 
for i=1:length(stims) 
    stimep = pop_epoch(ica,{stims{i}},[-.3,1]) ; 
    meaneps(i,:,:) = squeeze(mean(stimep.data(postelecs,:,:),1)) ; 
    stimep.data = stimep.data(:,:,1:90) ; 
    reseps = reshape(stimep.data,[64,333*90]) ; 
    for f=1:length(freqs)
        filt = eegfiltfft(reseps,ica.srate,freqs(f)-5,freqs(f)+5) ; 
        alldats(i,:,f,:,:) = reshape(filt,size(stimep.data)) ; 
    end
end
alldats = abs(alldats) ; 
stdtrials = squeeze(mean(std(mean(alldats,3),0,4),2)) ; 
[sv,si] = sort(stdtrials,2,'ascend') ; 
ntrials = 80 ; 
mdats = zeros(size(alldats,1),size(alldats,2),size(alldats,3),ntrials) ; 
for i=1:28
    mdats(i,:,:,:) = squeeze(mean(alldats(i,:,:,stimep.times<1000 & stimep.times>0,si(i,1:ntrials)),4) ...
                     - mean(alldats(i,:,:,stimep.times<0 & stimep.times>-250,si(i,1:ntrials)),4)) ;
end
figure,
for i=1:28 ; 
    subplot(4,7,i)
    topoplot(squeeze(mean(mean(mdats(i,:,6,:),3),4)),ica.chanlocs,'maplimits',[-.15,.15]) ;  title(i) ; 
end

gammadats = squeeze(mean(mean(mdats(:,:,5:7,:),3),4)) ; clim = [-.2,.2] ; 
subplot(2,3,1) ; topoplot((gammadats(1,:)),ica.chanlocs,'maplimits',clim) ; 
subplot(2,3,2) ; topoplot((gammadats(2,:)+gammadats(3,:)),ica.chanlocs,'maplimits',clim) ; 
subplot(2,3,3) ; topoplot((gammadats(4,:)+gammadats(5,:)),ica.chanlocs,'maplimits',clim) ; 
subplot(2,3,4) ; topoplot(sum(gammadats(6:13,:),1),ica.chanlocs,'maplimits',clim) ; 
subplot(2,3,5) ; topoplot(sum(gammadats(22:25,:),1),ica.chanlocs,'maplimits',clim) ; 
subplot(2,3,6) ; topoplot(sum(gammadats(14:21,:),1),ica.chanlocs) ; 









