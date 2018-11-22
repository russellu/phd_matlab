clear all ; close all; 
cd E:\angelina_2; 
eegs = dir('*gamma*vhdr'); 
for i=1:length(eegs)
eeg = pop_loadbv('.',eegs(i).name) ; 
eeg = pop_chanedit(eeg,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
eeg = pop_resample(eeg,256); 
if i==1; merged = eeg ; else merged=  pop_mergeset(merged,eeg); end
end

merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61); 
filtmerged = eegfiltfft(merged.data,merged.srate,1,90); 
trigs = {'S  2'}; 
%mergeica = merged; mergeica = pop_epoch(mergeica,trigs,[-1,10]); 
%mergeica = pop_runica(mergeica,'runica','maxsteps',128); 
[weights,sphere]= runica(filtmerged,'maxsteps',128); 
winv = pinv(weights*sphere); 
newmerged = merged; newmerged.data = weights*sphere*merged.data; 
clear ersp
for i=1:length(trigs);  disp(i); 
    epica = pop_epoch(newmerged,{trigs{i}},[-1,9]); 
    for j=1:64
            [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(j,:,:)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',200) ; 
    end
end
for i=1:64 ; subplot(5,13,i),imagesc(squeeze(mean(ersp(:,i,:,:),1)),[-3,3]); title(i);  end

% printshop@ubishops.ca
clear bersp
for i=1:length(trigs);  disp(i); 
    epica = pop_epoch(newmerged,{trigs{i}},[-1,9]); 
    for j=1:64
        for k=1:180
            [bersp(i,j,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(j,:,k)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off','timesout',200) ; 
        end
    end
end
mbersp = bersp - repmat(mean(bersp(1,:,:,:,times<0),5),[1,1,1,1,200]); 
cbersp = squeeze(mean(mbersp(:,[4,5,48],:,:,:),2)); 

erspmag = squeeze(mean(cbersp(31:120,:,:),1));
erspnomag = squeeze(mean(cbersp([1:30,121:end],:,:),1)); 

subplot(2,2,1); imagesc(times,freqs,erspnomag,[-3,3]); axis xy; vline([0,8],'k') ; xlabel('time(s)'); ylabel('frequency(hz)'); hline([50,70],'r'); title('no magnet');
subplot(2,2,2); imagesc(times,freqs,erspmag,[-3,3]); axis xy; vline([0,8],'k') ; xlabel('time(s)'); ylabel('frequency(hz)'); hline([50,70],'r'); title('magnet'); 

magtrials = squeeze(mean(mean(cbersp(31:120,25:35,times>0 & times<8),2),3)); 
nomagtrials = squeeze(mean(mean(cbersp([1:30,121:end],25:35,times>0 & times<8),2),3)); 
[h,p,ci,stats] = ttest(nomagtrials,magtrials); 
alltrials = [nomagtrials,magtrials]; 
subplot(2,2,3); barwitherr(std(alltrials,0,1)./sqrt(90),mean(alltrials,1)); set(gca,'XTickLabel',{'no magnet','magnet'}); 
title(['GAMMA differences: t = ',num2str(stats.tstat),', p=',num2str(p)]); 

magtrials = squeeze(mean(mean(cbersp(31:120,5:10,times>0 & times<8),2),3)); 
nomagtrials = squeeze(mean(mean(cbersp([1:30,121:end],5:10,times>0 & times<8),2),3)); 
[h,p,ci,stats] = ttest(nomagtrials,magtrials); 
alltrials = [nomagtrials,magtrials]; 
subplot(2,2,4); barwitherr(std(alltrials,0,1)./sqrt(90),mean(alltrials,1)); set(gca,'XTickLabel',{'no magnet','magnet'}); 
title(['ALPHA/BETA differences: t = ',num2str(stats.tstat),', p=',num2str(p)]); 





