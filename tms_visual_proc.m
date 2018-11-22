cd E:\data_for_Russ
savenames = {'ec_01.set','eo_01.set'};
for st=1:length(savenames)
eeg  = pop_loadset(savenames{st});
eeg = pop_chanedit(eeg,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;

tmsevents = find(strcmpi('TMS',{eeg.urevent.type}));
lats = cell2mat({eeg.urevent.latency}); 
tmslats = lats(tmsevents);
nepochs{st} = length(tmslats); 
for i=1:length(tmslats)
   eeg.data(:,tmslats(i)-25:tmslats(i)+25) = repmat(eeg.data(:,tmslats(i)-26)/2+eeg.data(:,tmslats(i)+26)/2,[1,51]);  
end


%figure,plot(eeg.data(46,:)); 

if st==1; merged = eeg;  else merged = pop_mergeset(merged,eeg); end
end
merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61); 

filtmerged = eegfiltfft(merged.data,eeg.srate,1,125); 
[weights,sphere] = runica(filtmerged,'maxsteps',128); 
winv = pinv(weights*sphere); 
for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),eeg.chanlocs); title(i); end

lfilt = eegfiltfft(merged.data,eeg.srate,0.2,125); 

acts = weights*sphere*lfilt; 
neweeg = merged;  neweeg.data = acts; 
allep = pop_epoch(neweeg,{'TMS'},[-.8,2]); 
clear mstersp; 
 for j=1:64 ; disp(j); 
    for k=1:size(allep.data,3)
        [mstersp(j,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data((j),:,k)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',80,'baseline',NaN,'verbose','off','timesout',200) ; 
    end
 end
bersp = mstersp - repmat(mean(mstersp(:,:,:,times<-.3),4),[1,1,1,200]); 
 
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(times,freqs,squeeze(mean(bersp(i,end-50:end,:,:),2)),[-3,3]); axis xy ; title(i); colormap jet; end

 
comps = [5,6,10]; 
 for i=1:3
    subplot(1,3,i); 
     plot(allep.times,squeeze(mean(allep.data(comps(i),:,1:50),3))); hold on; 
    plot(allep.times,squeeze(mean(allep.data(comps(i),:,end-50:end),3)),'r'); xlabel('time(ms)'); ylabel('component amp.'); title(['component ',num2str(comps(i))]);
 end
legend('EC','EO'); 
 
 
 
 
 
 