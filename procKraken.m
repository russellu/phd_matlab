% process kraken data
% run using eeglab 8.0 but will probably work on later versions
clear all ; close all ; 
cd c:/shared/kraken_Pilot/ ; ls 
EEG = pop_loadbv('.','M1.vhdr') ; 
EEG = pop_resample(EEG,250) ;
EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
%[chanpower,chanfreqs] = spectopo(EEG.data,0,EEG.srate,'plot','off') ; 
%figure,imagesc(chanpower) ; title('channel power') ;

% pre-process for ICA:
% here, i remove line noise and filter between 1 and 125Hz to improve the
% ICA decomposition (low frequencies mess it up). you can also play with
% the range, ie, if you are interested solely in mu, you may want to filter
% between 8-16Hz before running ICA, to make it more sensitive to that band.
filtEEG = EEG ; % copy the dataset, to preserve the unfiltered data for later use
filtEEG.data = filtEEG.data - eegfiltfft(filtEEG.data,filtEEG.srate,59.5,60.5) ; % filter out line noise
filtEEG.data = eegfiltfft(filtEEG.data,filtEEG.srate,1,125) ; % high pass filter above 1Hz

% triggers
trigs{1} = {'S  2','R  2','S  8'} ;
trigs{2} = {'S 32','R  8','S128'} ; 
trigs{3} = {'S 10','R 32','S 34'} ; 
trigs{4} = {'S130','R128','S 40'} ; 

ep = pop_epoch(filtEEG,[trigs{1},trigs{2},trigs{3},trigs{4}],[-2,2]) ; % all the epochs at once, for the decomposition
% below i run ICA on epoched data which i something i don't normally do,
% but there were alot of large transients in your data which made the
% continuous decomposition sub-optimal

% run ICA
ica = pop_runica(ep,'runica') ; 

% visualize the components => to see if ICA did a good job, you should see
% most of the components 1-20 looking like dipoles or monopoles. if you
% don't, something is wrong either you have a bad channel, or you have a
% lot of large transients in your data, or you forgot to high pass the data
% first.
for i=1:64 ; subplot(5,13,i) ; 
    topoplot(squeeze(ica.icawinv(:,i)),ica.chanlocs,'electrodes','on') ; title(i) ; 
end

%%% apply the weights to the continuous data, so you can have all
%%% frequencies in your components, and later re-epoch the ICA data as you wish
% there might be a function for this that i'm missing because its like 10
% lines for a simple task
EEG.icaact = icaact(EEG.data,ica.icaweights*ica.icasphere,0) ; 
EEG.icachansind = ica.icachansind ; 
EEG.icasphere = ica.icasphere ; 
EEG.icasplinefile = ica.icasplinefile ; 
EEG.icaweights = ica.icaweights ; 
EEG.icawinv = ica.icawinv ; 
EEG.data = EEG.data - eegfiltfft(EEG.data,ica.srate,59.5,60.5) ; % remove 60Hz in channel 
EEG.icaact = EEG.icaact - eegfiltfft(EEG.icaact,ica.srate,59.5,60.5) ; % remove 60Hz in component

% do some analysis now that all the preprocessing is done:
% first thing to check: ERSP in all conditions, baseline corrected 1s
% before go. first i compute the raw power in each single trial
% ('baseline',NaN), and after i compute the raw power i baseline correct
% and put the baseline corrected trials in "bersp"

% loop over all conditions:
for i=1:length(trigs)
    epochi = pop_epoch(EEG,trigs{i}(1),[-1,6]) ; % 7 second epoch according to go (index 1 in the trig vector for condition i)
    % computer ERSP over all components and trials:
    for comp=1:size(epochi.icaact,1) 
        for trial=1:size(epochi.icaact,3) ; 
            [ersp(i,comp,trial,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochi.icaact(comp,:,trial)),epochi.pnts,[epochi.xmin,epochi.xmax],epochi.srate,0,...
                'plotersp','off','plotitc','off','baseline',NaN,'winsize',64,'freqs',[1,120],'nfreqs',60,'timesout',200) ;          
        end
    end
end
bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,size(ersp,5)]) ; % baseline correction

% visualize all the ERSP from all conditions: a new figure for each
% condition, a new subplot for each component
figure,
for i=1:size(bersp,1) 
    figure,
    for j=1:64 ; subplot(5,13,j) ; 
        imagesc(times,freqs,squeeze(mean(bersp(i,j,:,:,:),3)),[-10,10]) ; axis xy ; title(['comp=',num2str(j)])
    end
end

%%%% now we do stats. for this type of experiment, i just do a simple
%%%% unpaired t-test between band and noband.


% first i resize all the ERSP so i can feed it into matlabs 'ttest2' as a
% single vector for each condition
clear conds
for i=1:size(bersp,1) 
    for j=1:size(bersp,3)
        conds(i,j,:) = reshape(squeeze(bersp(i,:,j,:,:)),[1,size(bersp,2)*size(bersp,4)*size(bersp,5)]) ; 
        
    end
end
% perform the t-test
[h1,p1,ci1,stats1] = ttest2(squeeze(conds(2,:,:)),squeeze(conds(1,:,:))) ; 
[h2,p2,ci2,stats2] = ttest2(squeeze(conds(4,:,:)),squeeze(conds(3,:,:))) ; 

% reshape the t-tests back to their original shape
tconds(1,:,:,:) = reshape(stats1.tstat,[size(bersp,2),size(bersp,4),size(bersp,5)]) ; 
tconds(2,:,:,:) = reshape(stats2.tstat,[size(bersp,2),size(bersp,4),size(bersp,5)]) ; 

pvals(1,:,:,:) = reshape(p1,[size(bersp,2),size(bersp,4),size(bersp,5)]) ; 
pvals(2,:,:,:) = reshape(p2,[size(bersp,2),size(bersp,4),size(bersp,5)]) ; 

% plot the t-value spectrograms
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(times,freqs,squeeze(tconds(2,i,:,:)),[-8,8]) ; axis xy ; title(['comp=',num2str(i)]) ; end


cbar = [-7,7] ; 
subplot(2,2,1) ; imagesc(times,freqs,squeeze(mean(bersp(3,3,:,:,:),3)),cbar) ; axis xy ; vline(0,'k') ; xlabel('time(s)') ; ylabel('freq(hz)' ); title('noband') ; 
subplot(2,2,2) ; imagesc(times,freqs,squeeze(mean(bersp(3,12,:,:,:),3)),cbar) ; axis xy ; vline(0,'k') ;xlabel('time(s)') ; ylabel('freq(hz)' ); title('noband') ; 
subplot(2,2,3) ; imagesc(times,freqs,squeeze(mean(bersp(4,3,:,:,:),3)),cbar) ; axis xy ; vline(0,'k') ; xlabel('time(s)') ; ylabel('freq(hz)' ); title('band') ; 
subplot(2,2,4) ; imagesc(times,freqs,squeeze(mean(bersp(4,12,:,:,:),3)),cbar) ; axis xy ; vline(0,'k') ;xlabel('time(s)') ; ylabel('freq(hz)' ); title('band') ; 

subplot(2,1,1) ; topoplot(squeeze(ica.icawinv(:,3)),ica.chanlocs,'electrodes','on') ; 
subplot(2,1,2) ; topoplot(squeeze(ica.icawinv(:,12)),ica.chanlocs,'electrodes','on') ; 




comps = [12,12] ; 

subplot(1,2,1) ; topoplot(squeeze(ica.icawinv(:,comps(1))),ica.chanlocs,'electrodes','labels') ; title('component 1 weight map') ;
subplot(1,2,2) ; topoplot(squeeze(ica.icawinv(:,comps(2))),ica.chanlocs,'electrodes','labels') ; title('component 2 weight map') ; 

startind = find(abs(times)==min(abs(times))) ; 

subplot(2,2,1) ; 
errorbar(squeeze(mean(mean(mean(bersp(1,comps,:,freqs>9 & freqs<14,:),2),3),4)),squeeze(std(mean(mean(bersp(1,comps,:,freqs>9 & freqs<14,:),2),4),0,3))./sqrt(55),'b') ; hold on ; 
errorbar(squeeze(mean(mean(mean(bersp(2,comps,:,freqs>9 & freqs<14,:),2),3),4)),squeeze(std(mean(mean(bersp(2,comps,:,freqs>9 & freqs<14,:),2),4),0,3))./sqrt(55),'r') ; 
xlim([1,200]) ; set(gca,'XTick',1:13:length(times),'XTickLabel',round(times(1:13:end)*1000)) ; xlabel('time(ms)') ; ylabel('9-14Hz power(db)') ; vline(startind,'g') ; hline(0,'k') ; title('mu ERD, 2s') ; 
legend({'noband','band'}) ; 
subplot(2,2,2) ; imagesc(times,freqs,squeeze(mean(tconds(1,comps,:,:),2)),[-10,10]) ; axis xy ; vline(0,'k') ; colorbar ; title('t values band-noband (2s)') ; xlabel('time(s)') ; ylabel('freq(hz)') ; 

subplot(2,2,3) ; 
errorbar(squeeze(mean(mean(mean(bersp(3,comps,:,freqs>9 & freqs<14,:),2),3),4)),squeeze(std(mean(mean(bersp(3,comps,:,freqs>9 & freqs<14,:),2),4),0,3))./sqrt(55),'b') ; hold on ; 
errorbar(squeeze(mean(mean(mean(bersp(4,comps,:,freqs>9 & freqs<14,:),2),3),4)),squeeze(std(mean(mean(bersp(4,comps,:,freqs>9 & freqs<14,:),2),4),0,3))./sqrt(55),'r') ; 
xlim([1,200]) ; set(gca,'XTick',1:13:length(times),'XTickLabel',round(times(1:13:end)*1000)) ; xlabel('time(ms)') ; ylabel('9-14Hz power(db)') ; vline(startind,'g') ; hline(0,'k') ; title('mu ERD, 4s') ; 
legend({'noband','band'}) ; 
subplot(2,2,4) ; imagesc(times,freqs,squeeze(mean(tconds(2,comps,:,:),2)),[-10,10]) ; axis xy ; vline(0,'k') ; colorbar ; title('t values band-noband (4s)') ; xlabel('time(s)') ; ylabel('freq(hz)') ; 





plot(squeeze(mean(mean(bersp(3,[3,14],:,freqs>10 & freqs<14,:),2),4))','r'); hold on ; 

plot(squeeze(mean(mean(bersp(4,[3,14],:,freqs>10 & freqs<14,:),2),4))','b');
ep2 = pop_epoch(EEG,[trigs{1},trigs{2},trigs{3},trigs{4}],[-4,0]) ;
[s,f] = spectopo(ep2.icaact,0,EEG.srate,'plot','off') ;
plot(mean(s([12],1:length(f)/2),1),'b','LineWidth',2); set(gca,'XTick',1:25:length(f)/2,'XTickLabel',round(f(1:25:length(f)/2))) ; xlabel('frequency(hz)') ; ylabel('power(log)') ; hold on ; 
plot(mean(s([14],1:length(f)/2),1),'r','LineWidth',2); set(gca,'XTick',1:25:length(f)/2,'XTickLabel',round(f(1:25:length(f)/2))) ; xlabel('frequency(hz)') ; ylabel('power(log)') ; 
legend({'contralateral baseline power','ipsilateral baseline power'}) ; 


