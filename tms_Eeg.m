clear all ; close all; 

cd('E:\TMS EEG');
eeg = pop_loadbv('.','C01_SC_TMS_TS.vhdr');
eeg = pop_chanedit(eeg,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;

stimtrigs = {'S  1','S  2','S  4','S  8'};
%stimtrigs = {'S  4','S  5','S  6','S 12'};
neweeg = eeg; 
for tr=1:length(stimtrigs)
    lats = cell2mat({eeg.urevent.latency});
    trigs = {eeg.urevent.type};
    trig_inds = find(strcmpi(stimtrigs{tr},trigs)); 
    clear trigdat triginds; 
    for i=1:length(trig_inds)
       trigdat(i,:,:) = eeg.data(:,lats(trig_inds(i))-12:lats(trig_inds(i))+12);  
       triginds(i,:) = lats(trig_inds(i))-12:lats(trig_inds(i))+12; 
    end
    
    for i=1:64 
        corrs = corr(squeeze(trigdat(:,i,:))'); 
        [sv,si] = sort(corrs,2,'descend'); 
        for j=1:length(si)
            neweeg.data(i,triginds(j,:)) = eeg.data(i,triginds(j,:)) - squeeze(mean(trigdat(si(j,1),i,:),1))';          
        end
    end
    
end

reseeg = pop_resample(neweeg,250);
resfilt = reseeg.data - eegfiltfft(reseeg.data,reseeg.srate,0,1) - eegfiltfft(reseeg.data,reseeg.srate,59.5,60.5); 
[weights,sphere] = runica(resfilt,'maxsteps',128);
winv = pinv(weights*sphere); 
figure,for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),eeg.chanlocs); title(i); end
acts = weights*sphere*resfilt; 

clear allersp; 
icaeeg = reseeg; icaeeg.data = acts; 
for s=1%:length(stimtrigs) ; disp(s) ; 
    allep = pop_epoch(icaeeg,{stimtrigs{s}},[-1,3]) ; 
for i=1:64 
        [allersp(s,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',130,'baseline',0,'verbose','off','timesout',200) ; 
end
end   

ep1 = pop_epoch(icaeeg,{stimtrigs{1}},[-1,2]) ; 
ep2 = pop_epoch(icaeeg,{stimtrigs{2}},[-1,2]) ; 
ep3 = pop_epoch(icaeeg,{stimtrigs{3}},[-1,2]) ; 
ep4 = pop_epoch(icaeeg,{stimtrigs{4}},[-1,2]) ; 

colors = {'r','b','g','m'};

figure,tinds = find(ep1.times>-100 & ep1.times<500); 
subplot(3,2,1); 
topoplot(winv(:,13),eeg.chanlocs); title('component 16'); 
subplot(3,2,2); 
plot(ep1.times(tinds),squeeze(mean(ep1.data(13,tinds,:),3))); hold on ;
plot(ep1.times(tinds),squeeze(mean(ep2.data(13,tinds,:),3)));
plot(ep1.times(tinds),squeeze(mean(ep3.data(13,tinds,:),3)));
plot(ep1.times(tinds),squeeze(mean(ep4.data(13,tinds,:),3))); vline(0,'k'); xlabel('time(ms)') ;ylabel('component amplitude'); 
legend({'S  1','S  2','S  4','S  8'});title('component 16'); 

subplot(3,2,3); 
topoplot(winv(:,16),eeg.chanlocs); title('component 13');
subplot(3,2,4); 
plot(ep1.times(tinds),squeeze(mean(ep1.data(16,tinds,:),3))); hold on ;
plot(ep1.times(tinds),squeeze(mean(ep2.data(16,tinds,:),3)));
plot(ep1.times(tinds),squeeze(mean(ep3.data(16,tinds,:),3)));
plot(ep1.times(tinds),squeeze(mean(ep4.data(16,tinds,:),3))); vline(0,'k'); xlabel('time(ms)') ;ylabel('component amplitude'); 
legend({'S  1','S  2','S  4','S  8'});title('component 13');

subplot(3,2,5); 
topoplot(winv(:,7),eeg.chanlocs); title('component 7');
subplot(3,2,6); 
plot(ep1.times(tinds),squeeze(mean(ep1.data(7,tinds,:),3))); hold on ;
plot(ep1.times(tinds),squeeze(mean(ep2.data(7,tinds,:),3)));
plot(ep1.times(tinds),squeeze(mean(ep3.data(7,tinds,:),3)));
plot(ep1.times(tinds),squeeze(mean(ep4.data(7,tinds,:),3))); vline(0,'k'); xlabel('time(ms)') ;ylabel('component amplitude'); 
legend({'S  1','S  2','S  4','S  8'});title('component 7');
%{
ts = 311409-10:311409+10; 
for i=1:length(ts) ; subplot(4,5,i) ; topoplot(eeg.data(:,ts(i)),eeg.chanlocs,'maplimits',[0,4000]); end

%trigs = unique({eeg.urevent.type});
trigs = {'S  1','S  2','S  4'};
filteeg = eeg.data - eegfiltfft(eeg.data,eeg.srate,0,1) - eegfiltfft(eeg.data,eeg.srate,59.5,60.5);
[weights,sphere] = runica(filteeg,'maxsteps',128);
winv = pinv(weights*sphere); 
for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),eeg.chanlocs); title(i); end
acts = weights*sphere*filteeg; 

for tr=1:length(trigs)
    eps = pop_epoch(eeg,{trigs{tr}},[-1,1]); 
    subplot(1,3,tr) ; plot(squeeze(mean(eps.data,3))') ;
end



ep1 = pop_epoch(eeg,{stimtrigs{1}},[-.1,.1]) ; 
ep2 = pop_epoch(eeg,{stimtrigs{2}},[-.1,.1]) ; 
ep3 = pop_epoch(eeg,{stimtrigs{3}},[-.1,.1]) ; 
ep4 = pop_epoch(eeg,{stimtrigs{4}},[-.1,.1]) ; 

plot(ep1.times,squeeze(mean(mean(ep1.data,1),3))); hold on; 
plot(ep1.times,squeeze(mean(mean(ep2.data,1),3))); 
plot(ep1.times,squeeze(mean(mean(ep3.data,1),3))); 
plot(ep1.times,squeeze(mean(mean(ep4.data,1),3))); xlabel('time(ms)') ;ylabel('\muV'); 
legend({'S  1','S  2','S  4','S  8'});title('TMS pulses');


%}






clear allersp; 
neweeg = reseeg; neweeg.data = acts; 
for s=1%:length(stimtrigs) ; disp(s) ; 
    allep = pop_epoch(neweeg,{stimtrigs{s}},[-1,1.5]) ; 
    allep.data(:,allep.times>-25 & allep.times<25,:) = 0; 
for i=1:64 
        [allersp(s,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',100) ; 
end
end   












