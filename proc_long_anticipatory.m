clear all ; close all; 
cd c:/shared/long_anticipatory; 
eegs = dir('*vhdr'); 
for i=1:length(eegs)
eeg = pop_loadbv('.',eegs(i).name) ; 
eeg = pop_chanedit(eeg,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
eeg = pop_resample(eeg,256); 
if i==1; merged = eeg ; else merged=  pop_mergeset(merged,eeg); end
end

merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61); 
merged.data(32,:) = rand(1,size(merged.data,2))/10; 
filtmerged = eegfiltfft(merged.data,merged.srate,1,90); 
trigs = {'S 11','S 12','S 14'}; 
%mergeica = merged; mergeica = pop_epoch(mergeica,trigs,[-1,10]); 
%mergeica = pop_runica(mergeica,'runica','maxsteps',128); 
[weights,sphere]= runica(filtmerged,'maxsteps',128); 
winv = pinv(weights*sphere); 
newmerged = merged; newmerged.data = weights*sphere*merged.data; 
clear ersp
for i=1:length(trigs);  disp(i); 
    epica = pop_epoch(newmerged,{trigs{i}},[-5,5]); 
    for j=1:64
        for t=1:size(epica.data,3)
            [ersp(i,j,t,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(j,:,t)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off','timesout',200) ; 
        end
    end
end
bersp = ersp - repmat(mean(ersp(:,:,:,:,times<-3),5),[1,1,1,1,200]); 
for i=1:64 ; subplot(5,13,i),imagesc(squeeze(mean(mean(bersp(3,i,:,:,:),1),3)),[-3,3]); title(i); axis xy; end

comps = [5,19]; 
trials = {1:40,41:80,81:100,101:120}; titles = {'long ISI','short ISI','0.5s focus','3s focus'};
for i=1:4 ; subplot(3,4,i) ; imagesc(times,freqs,squeeze(mean(mean(bersp(1,comps,trials{i},:,:),2),3)),[-5,5]); 
    axis xy; title(titles{i});  
    if i==1 || i==2
        vline([0,-1.5],'k'); 
    elseif i==3
        vline([0,-0.5],'k'); 
    elseif i==4
        vline([0,-3],'k'); 
    end
end
for i=1:4 ; subplot(3,4,i+4) ; imagesc(times,freqs,squeeze(mean(mean(bersp(2,comps,trials{i},:,:),2),3)),[-5,5]); 
    axis xy; title(titles{i}); 
        if i==1 || i==2
        vline([0,-1.5],'k'); 
    elseif i==3
        vline([0,-0.5],'k'); 
    elseif i==4
        vline([0,-3],'k'); 
    end
end
for i=1:4 ; subplot(3,4,i+8) ; imagesc(times,freqs,squeeze(mean(mean(bersp(3,comps,trials{i},:,:),2),3)),[-5,5]); 
    axis xy; title(titles{i}); 
        if i==1 || i==2
        vline([0,-1.5],'k'); 
    elseif i==3
        vline([0,-0.5],'k'); 
    elseif i==4
        vline([0,-3],'k'); 
    end
end


