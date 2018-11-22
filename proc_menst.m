clear all ; close all; 
cd E:\angelina_2; 
eegs = dir('*_luteal_gamma*vhdr'); 
for i=1:length(eegs)
eeg = pop_loadbv('.',eegs(i).name) ; 
eeg = pop_chanedit(eeg,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
eeg = pop_resample(eeg,256); 
if i==1; merged = eeg ; else merged=  pop_mergeset(merged,eeg); end
end

merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61); 
filtmerged = eegfiltfft(merged.data,merged.srate,1,90); 
trigs = {'S  1','S  2'}; 

[weights,sphere]= runica(filtmerged,'maxsteps',128); 
winv = pinv(weights*sphere); 
newmerged = merged; newmerged.data = weights*sphere*merged.data; 
clear ersp
for i=1:length(trigs);  disp(i); 
    epica = pop_epoch(newmerged,{trigs{i}},[-2,17]); 
    for j=1:64
            [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(j,:,:)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',200) ; 
    end
end
for i=1:64 ; subplot(5,13,i),imagesc(squeeze(mean(ersp(:,i,:,:),1)),[-8,8]); axis xy ; colormap jet; title(i);  end
comps = [2,4,7,22]; 

subplot(1,2,1) ; imagesc(times,freqs,squeeze(mean(ersp(1,comps,:,:),2)),[-5,5]); axis xy; colormap jet; title('high contrast'); 
subplot(1,2,2) ; imagesc(times,freqs,squeeze(mean(ersp(2,comps,:,:),2)),[-5,5]); axis xy; colormap jet; title('low contrast'); xlabel('time(s)'); ylabel('frequency(hz)'); 




