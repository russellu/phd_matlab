clear all ; close all ; 
cd C:\shared\russ_tuning ;   
discr = dir('*allstim*vhdr') ; 
for i=1:length(discr)
   EEG = pop_loadbv('.',discr(i).name) ; 
   EEG = pop_resample(EEG,256) ; 
   if i==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG) ; end
end
merged = pop_chanedit(merged,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
temp = merged.data(1:32,:) ; merged.data(1:32,:) = merged.data(33:64,:) ; merged.data(33:64,:) = temp ; 
merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61) - eegfiltfft(merged.data,merged.srate,84,86); 
stims = {'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8','S  9','S 10','S 11','S 12','S 13','S 14','S 15','S 16','S 17','S 18','S 19','S 20','S 21','S 22',...
    'S 23','S 24','S 25','S 26','S 27','S 28','S 41','S 42','S 43','S 44','S 45','S 46','S 47','S 48','S 49','S 50','S 51','S 52','S 53','S 54','S 55','S 56'} ;
mergefilt = merged ; mergefilt.data = eegfiltfft(merged.data,merged.srate,40,70) ; %mergep = pop_epoch(mergefilt,stims,[-.3,1.3]) ; 
%mergedat = reshape(mergep.data,[64,410*4400]) ; [weights,sphere] = runica(mergefilt.data(:,1:5:end),'maxsteps',128) ; winv = pinv(weights*sphere) ; 
ica = pop_runica(mergefilt,'runica') ; ica.icaact = ica.icaweights*ica.icasphere*merged.data ; 
ica.data = ica.icaact ; 

clear allersp ; 
for s=1:length(stims) ; disp(s) ; 
    allep = pop_epoch(ica,{stims{s}},[-.3,1.3]) ;
    for i=1:64 ;
            [allersp(s,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off') ; 
                
    end
end
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(allersp(:,i,:,:),1)),[-3,3]) ; title(i) ; end
goodcs = [15,27] ; 

mc = squeeze(mean(mean(allersp(:,goodcs,:,times>.25 & times<1),2),4)) ; 
subplot(2,2,1) ; low = squeeze(mean(mc(29:end,freqs>8 & freqs<25),2)) ; high = squeeze(mean(mc(29:end,freqs>40 & freqs<80),2)) ; 
plot(low,high,'o') ; lsline ; title(num2str(corr2(low,high))) ; %title('orientation tuning correlation (alpha/beta vs gamma)') ; 
subplot(2,2,2) ; low = squeeze(mean(mc(1:28,freqs>8 & freqs<25),2)) ; high = squeeze(mean(mc(1:28,freqs>40 & freqs<80),2)) ; 
plot(low,high,'o') ; lsline ; title(num2str(corr2(low,high))) ; %title('retinotopic tuning correlation (alpha/beta vs gamma)') ; 


clear tersp
for s=1:length(stims) ; disp(s) ; 
    allep = pop_epoch(ica,{stims{s}},[-.25,1.25]) ;
    for i=1:length(goodcs) ;
        for j=1:70
            [tersp(s,i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(goodcs(i),:,j)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off','timesout',100) ; 
        end       
    end
end
btersp = tersp - repmat(mean(tersp(:,:,:,:,times<0),5),[1,1,1,1,100]) ; 
mbtersp = squeeze(mean(mean(btersp(:,:,:,:,times>.5 & times<1),2),5)) ; 

errorbar(squeeze(mean(mean(mbtersp(:,:,freqs>40 & freqs<60),2),3)),squeeze(std(mean(mbtersp(:,:,freqs>40 & freqs<60),3),0,2))/sqrt(100)) ; 
errorbar(squeeze(mean(mean(mbtersp(:,:,freqs>8 & freqs<25),2),3)),squeeze(std(mean(mbtersp(:,:,freqs>8 & freqs<25),3),0,2))/sqrt(100)) ; 
f60 = 1:60 ; 
clear allts ; 
for i=1:38
    for j=1:38
        for f=1:60 ; 
            [h,p,ci,stats] = ttest(squeeze(mean(mbtersp(i,:,(f60>f-3 & f60<f+3)),3)),squeeze(mean(mbtersp(j,:,(f60>f-3 & f60<f+3)),3))) ; 
            allts(i,j,f) = stats.tstat ; 
        end 
    end
end
allts(isnan(allts)) = 0 ; 

