clear all ; close all; 
cd c:/shared/greg_pupil_eeg ; 
eegs = dir('*vhdr'); 
for i=1:length(eegs)
eeg = pop_loadbv('.',eegs(i).name) ; 
eeg = pop_chanedit(eeg,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
eeg = pop_resample(eeg,256); 
if i==1; merged = eeg ; else merged=  pop_mergeset(eeg,merged); end
end

merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61); 
merged.data(32,:) = rand(1,size(merged.data,2))/10; 
filtmerged = eegfiltfft(merged.data,merged.srate,1,128); 
trigs = {'S 11','S 12','S 13','S 14'}; 

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


clear bersp
for i=1:length(trigs);  disp(i); 
    epica = pop_epoch(newmerged,{trigs{i}},[-2,3.5]); 
    alleps{i} = epica; 
    for j=1:64
        for t=1:80
            [bersp(i,j,t,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(j,:,t)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off','timesout',200) ; 
        end
    end
end
mbersp = bersp - repmat(mean(bersp(:,:,:,:,times>-.5 & times<0),5),[1,1,1,1,200]); 
comps = [5,19]; titles = {'unperturbed','33%contrast','60%rnd','no stim'};
for i=1:4 ; subplot(2,2,i) ;
    imagesc(times,freqs,squeeze(mean(mean(mbersp(i,comps,41:end,:,:),2),3)),[-6,6]); axis xy ; title(titles{i}); 
    vline(-1.5,'k'); vline(0,'k'); 
end

for i=1:4 ; subplot(2,2,i) ; plot(squeeze(mean(mean(alleps{i}.data(comps,:,41:end),1),3))); ylim([-15,-0]); title(titles{i});  end








