clear all ; close all; 
cd c:/shared/shortstim; 
eegs = dir('*vhdr'); 
for i=1:2%length(eegs)
eeg = pop_loadbv('.',eegs(i).name) ; 
eeg = pop_chanedit(eeg,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
eeg = pop_resample(eeg,256); 
if i==1; merged = eeg ; else merged=  pop_mergeset(eeg,merged); end
end

merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61); 
merged.data(32,:) = rand(1,size(merged.data,2))/10; 
filtmerged = eegfiltfft(merged.data,merged.srate,40,90); 
trigs = {'S 11','S 12','S 13','S 14','S 15','S 21','S 22','S 23','S 24','S 25','S 31','S 32','S 33','S 34','S 35'}; 
%mergeica = merged; mergeica = pop_epoch(mergeica,trigs,[-1,10]); 
%mergeica = pop_runica(mergeica,'runica','maxsteps',128); 
[weights,sphere]= runica(filtmerged,'maxsteps',128); 
winv = pinv(weights*sphere); 
newmerged = merged; newmerged.data = weights*sphere*merged.data; 
clear ersp
for i=1:length(trigs);  disp(i); 
    epica = pop_epoch(newmerged,{trigs{i}},[-.5,3.5]); 
    for j=1:64
            [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(j,:,:)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',200) ; 
    end
end
for i=1:64 ; subplot(5,13,i),imagesc(squeeze(mean(ersp(:,i,:,:),1)),[-3,3]); title(i);  end
comps = [10,16];
contlvls = [.05,.15,.25,.5,1]; 
rndlvls = [0.05,0.1,0.25,0.5,1]; 
circrads = [25,50,100,150,200]; 
subplot(2,2,1); 
plot(freqs(12:45),squeeze(mean(mean(ersp(1:5,comps,12:45,times>0.5 & times<3),2),4))','LineWidth',4);  
title('contrast'); xlabel('frequency(hz)') ; ylabel('power(db)'); hline(0,'k'); ylim([-1,6]); 
legend({num2str(contlvls(1)),num2str(contlvls(2)),num2str(contlvls(3)),num2str(contlvls(4)),num2str(contlvls(5))}); 
subplot(2,2,2); 
plot(freqs(12:45),squeeze(mean(mean(ersp(6:10,comps,12:45,times>0.5 & times<3),2),4))','LineWidth',4);  
title('randomization');ylabel('power(db)'); hline(0,'k'); ylim([-1,6]); 
legend({num2str(rndlvls(1)),num2str(rndlvls(2)),num2str(rndlvls(3)),num2str(rndlvls(4)),num2str(rndlvls(5))}); 
subplot(2,2,3); 
plot(freqs(12:45),squeeze(mean(mean(ersp(11:15,comps,12:45,times>0.5 & times<3),2),4))','LineWidth',4);  
title('size');ylabel('power(db)'); hline(0,'k'); ylim([-1,6]); 
legend({num2str(circrads(1)),num2str(circrads(2)),num2str(circrads(3)),num2str(circrads(4)),num2str(circrads(5))}); 

subplot(2,2,1);
bar(mersp(1:5)) ; title('contrast'); 
subplot(2,2,2);
bar(mersp(6:10));title('randomization'); 
subplot(2,2,3);
bar(mersp(11:15));title('size'); 

mersp = squeeze(mean(mean(ersp(:,comps,:,times>0.5 & times<3),2),4)); 



    
    