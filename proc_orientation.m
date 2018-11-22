clear all ; close all ; 
cd C:\shared\orientation
discr = dir('russ_orientation*vhdr') ; 
for i=1:length(discr)
   EEG = pop_loadbv('.',discr(i).name) ; 
   EEG = pop_resample(EEG,256) ; 
   if i==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG) ; end
end
merged = pop_chanedit(merged,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
temp = merged.data(1:32,:) ; merged.data(1:32,:) = merged.data(33:64,:) ; merged.data(33:64,:) = temp ; 
merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61) ; 
mergefilt = merged ; mergefilt.data = eegfiltfft(merged.data,merged.srate,40,65) ; 
ica = pop_runica(mergefilt,'runica') ; ica.icaact = ica.icaweights*ica.icasphere*merged.data ; 
ica.data = ica.icaact ; 
stims = {'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8','S  9','S 10','S 11','S 12','S 13','S 14','S 15','S 16','S 17','S 18','S 19','S 20','S 21','S 22','S 23','S 24'} ;
clear allersp ; 
for s=1:length(stims) ; disp(s) ; 
    allep = pop_epoch(ica,{stims{s}},[-.25,1.25]) ;
    for i=1:64 ;
            [allersp(s,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off') ; 
                
    end
end
for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(allersp(:,i,:,:),1)),[-8,8]) ; title(i) ; end
goodcs = [8,9] ; 
gamma = squeeze(mean(mean(mean(allersp(:,goodcs,freqs>50 & freqs<65,times>0 & times<1),2),3),4)) ; 
alpha = squeeze(mean(mean(mean(allersp(:,goodcs,freqs>8 & freqs<25,times>0 & times<1),2),3),4)) ;


clear tersp
for s=1:length(stims) ; disp(s) ; 
    allep = pop_epoch(ica,{stims{s}},[-.25,1.25]) ;
    for i=1:length(goodcs) ;
        for j=1:70
            [tersp(s,i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(goodcs(i),:,j)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off') ; 
        end       
    end
end

btersp = tersp - repmat(mean(tersp(:,:,:,:,times<0),5),[1,1,1,1,200]) ; 
mbtersp = squeeze(mean(mean(btersp(:,:,:,:,times>0.5 & times<1),2),5)) ; 
for i=1:24 ; subplot(4,8,i) ; imagesc(squeeze(mbtersp(i,:,:)),[-10,10]) ; end
for i=1:60 ; subplot(6,10,i) ,errorbar(squeeze(mean(mean(mbtersp(:,:,freqs>i-3 & freqs<i+3),2),3)),squeeze(std(mean(mbtersp(:,:,freqs>i-3 & freqs<i+3),3),0,2))./sqrt(70))
    xlim([0,24]) ; 
end
    angles = 1:15:360 ; 
subplot(2,2,2) ; 
shadedErrorBar([],squeeze(mean(mbtersp(1,:,:),2)),squeeze(std(mbtersp(1,:,:),0,2))./sqrt(70),{'b'}) ; hold on ; 
shadedErrorBar([],squeeze(mean(mbtersp(7,:,:),2)),squeeze(std(mbtersp(7,:,:),0,2))./sqrt(70),{'r'}) ; hline(0,'k'); 
set(gca,'XTick',0:5:60,'XTickLabel',0:10:120) ; xlabel('frequency(hz)') ; ylabel('power(db)') ; title('vertical(blue) vs horizontal(red)') ; 
subplot(2,2,1) ; imagesc(freqs,angles,squeeze(mean(mbtersp,2)),[-8,8]) ; xlabel('frequency(hz)') ; ylabel('orientation(deg)') ; 
subplot(2,2,3) ; errorbar(squeeze(mean(mean(mbtersp(:,:,freqs>40 & freqs<70),2),3)),squeeze(std(mean(mbtersp(:,:,freqs>40 & freqs<70),3),0,2))./sqrt(70))
set(gca,'XTick',1:3:24,'XTickLabel',angles(1:3:24)) ; title('gamma (50-70Hz) tuning') ; ylabel('power(db)') ; 
subplot(2,2,4) ; errorbar(squeeze(mean(mean(mbtersp(:,:,freqs>10 & freqs<25),2),3)),squeeze(std(mean(mbtersp(:,:,freqs>10 & freqs<25),3),0,2))./sqrt(70))
set(gca,'XTick',1:3:24,'XTickLabel',angles(1:3:24)) ; title('alpha/beta (8-25Hz) tuning') ; ylabel('power(db)') ; 
f60 = 1:60 ; 
clear allts ; 
for i=1:24
    for j=1:24
        for f=1:60 ; 
            [h,p,ci,stats] = ttest(squeeze(mean(mbtersp(i,:,find(f60>f-3 & f60<f+3)),3)),squeeze(mean(mbtersp(j,:,find(f60>f-3 & f60<f+3)),3))) ; 
            allts(i,j,f) = stats.tstat ; 
        end 
    end
end
allts(isnan(allts)) = 0 ; 
subplot(1,2,1) ; 
bar(squeeze(mean(mean(btersp(1,:,:,31,times>.5&times<1),2),5)))
subplot(1,2,2) ; 
bar(squeeze(mean(mean(btersp(7,:,:,31,times>.5&times<1),2),5)))

subplot(2,2,1) ; imagesc(freqs,freqs,corrs) ; title('correlation of tuning curves (all frequencies)') ; xlabel('frequency(hz)') ; ylabel('frequency(hz)');  
subplot(2,2,2) ; 
plot(squeeze(mean(mbtersp(:,:,16),2)),squeeze(mean(mbtersp(:,:,31),2)),'o') ; lsline ; xlabel('low gamma (32Hz)') ; ylabel('high gamma(62Hz)') ; 
title(num2str(corr2(squeeze(mean(mbtersp(:,:,16),2)),squeeze(mean(mbtersp(:,:,31),2)))))
subplot(2,2,3) ; errorbar(squeeze(mean(mbtersp(:,:,31),2)),squeeze(std(mbtersp(:,:,31),0,2))./sqrt(70),'b') ; hold on ; 
errorbar(squeeze(mean(mbtersp(:,:,16),2)),squeeze(std(mbtersp(:,:,16),0,2))./sqrt(70),'r') ; legend({'low gamma','high gamma'}) ; title('low/high gamma tuning curves') ;  
xlabel('orientation (A.U)') ; ylabel('power(db)') ; 




