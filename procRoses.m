% process kraken data
% run using eeglab 8.0 but will probably work on later versions
clear all ; close all ; 
cd c:/shared/FourRoses_Pilot/FourRoses_EEG_DATA ; ls 
EEG = pop_loadbv('.','M3.vhdr') ; 
EEG = pop_resample(EEG,250) ;
EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
%[chanpower,chanfreqs] = spectopo(EEG.data,0,EEG.srate,'plot','off') ; 
%figure,imagesc(chanpower) ; title('channel power') ; 

% pre-process for ICA:
filtEEG = EEG ; % copy the dataset, to preserve the unfiltered data for later use
filtEEG.data = filtEEG.data - eegfiltfft(filtEEG.data,filtEEG.srate,59.5,60.5) ; % filter out line noise
filtEEG.data = eegfiltfft(filtEEG.data,filtEEG.srate,1,120) ; % high pass filter above 1Hz

% get the triggers and use them to filter out bad epochs 
trigs{1} = {'S  2'} ;trigs{2} = {'S 32'} ; trigs{3} = {'R  8'} ; trigs{4} = {'S130'} ; trigs{5} = {'S136'} ; trigs{6} = {'R128'} ;
trigs{7} = {'S162'} ;trigs{8} = {'S170'} ; trigs{9} = {'R 34'} ; trigs{10} = {'S 40'} ; trigs{11} = {'S 42'} ; trigs{12} = {'R136'} ;

ep = pop_epoch(filtEEG,[trigs{1},trigs{2},trigs{3},trigs{4},trigs{5},trigs{6},trigs{7},trigs{8},trigs{9},trigs{10},trigs{11},trigs{12}],[-2,4]) ; % all the epochs at once, for the decomposition
% below i run ICA on epoched data which i something i don't normally do,
% but there were alot of large transients in your data which made the
% continuous decomposition sub-optimal

% run ICA
ica = pop_runica(ep,'runica') ; 

% visualize the components
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
% before go 


% loop over all conditions:
clear ersp 
for i=1:length(trigs)
    epochi = pop_epoch(EEG,trigs{i}(1),[-1,6]) ; % 7 second epoch according to go (index 1 in the trig vector for condition i)
    % computer ERSP over all components and trials:
    for comp=1:size(epochi.icaact,1) 
        for trial = 1:size(epochi.data,3)
            [ersp(i,comp,trial,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochi.icaact(comp,:,trial)),epochi.pnts,[epochi.xmin,epochi.xmax],epochi.srate,0,...
                'plotersp','off','plotitc','off','baseline',NaN,'winsize',64,'freqs',[1,120],'nfreqs',60,'timesout',100) ; 
            
        end     
    end
end

bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,100]) ; 


clear conds
for i=1:size(bersp,1) 
    for j=1:size(bersp,3)
        conds(i,j,:) = reshape(squeeze(bersp(i,:,j,:,:)),[1,size(bersp,2)*size(bersp,4)*size(bersp,5)]) ; 
        
    end
end

[h1,p1,ci1,stats1] = ttest2(squeeze(conds(10,:,:)),squeeze(conds(1,:,:))) ; 
tconds(1,:,:,:) = reshape(stats1.tstat,[size(bersp,2),size(bersp,4),size(bersp,5)]) ; 
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(tconds(1,i,:,:)),[-8,8]) ; title(i) ;  end

tstart = find(abs(times)==min(abs(times))) ; 

for i=1:384000 ; 
    [p,anovatab,stats] = anova1(squeeze(conds(:,:,i))',[],'off') ; disp(i/384000) ; 
    fs(i) = anovatab{2,5} ; 
end

resfs = reshape(fs,[64,60,100]) ; 
figure,for i=1:25 ; subplot(5,5,i) ; imagesc(times,freqs,squeeze(resfs(i,:,:)),[0,10]) ; title(['component ',num2str(i)]) ;axis xy ; end
figure,for i=1:25 ; subplot(5,5,i) ; topoplot(squeeze(ica.icawinv(:,i)),ica.chanlocs) ; title(['component ',num2str(i)]) ; end
figure,for i=1:25 ; subplot(5,5,i) ; imagesc(times,freqs,squeeze(mean(mean(bersp(:,i,:,:,:),1),3)),[-8,8]) ; title(['component ',num2str(i)]) ; axis xy ;  end

comps = [5,6,8] ; 

for c=1:25 ; 

fh = figure('Position',[0,0,1500,500]);subplot(1,3,1) ; 
topoplot(squeeze(ica.icawinv(:,c)),ica.chanlocs,'electrodes','on') ; 
flags = {'0','60','120','180','240','300','anti_0','anti_60','anti_120','anti_180','anti_240','anti_300'} ; 
subplot(1,3,2) ; 
plot(squeeze(mean(mean(bersp(1:6,c,:,freqs>8& freqs<16,:),3),4))','LineWidth',2) ; hline(0,'k') ; vline(tstart,'k') ;
legend(flags(1:6)) ; 
set(gca,'XTick',1:10:length(times),'XTickLabel',round(times(1:10:end)*1000)) ; ylabel('mu ERD (8-16Hz)') ; xlabel('time(ms)') ; 
subplot(1,3,3) ; 
plot(squeeze(mean(mean(bersp(7:12,c,:,freqs>8& freqs<16,:),3),4))','LineWidth',2) ; hline(0,'k') ; vline(tstart,'k') ; 
legend(flags(7:12)) ; 
set(gca,'XTick',1:10:length(times),'XTickLabel',round(times(1:10:end)*1000)) ; ylabel('mu ERD (8-16Hz)') ; xlabel('time(ms)') ; 

end


plot(squeeze(mean(tconds(1,6,freqs>8 & freqs<16,:),3)),'LineWidth',3) ; hline(0,'k') ; 
set(gca,'XTick',1:10:length(times),'XTickLabel',round(times(1:10:end)*1000)) ; ylabel('t-value') ; xlabel('time(ms)') ; 



for i=1:12 ; subplot(2,6,i) ; imagesc(times,freqs,squeeze(mean(bersp(i,12,:,:,:),3)),[-10,10]) ; axis xy ; vline(0,'k') ; xlabel('time(s)') ; ylabel('freq(hz)') ; end 


meanpre = squeeze(mean(bersp(:,:,:,:,times>3 & times<4.5),5)) ; 
subplot(1,2,1) ; barwitherr(squeeze(std(mean(meanpre(1:6,12,:,freqs>9 & freqs<15),4),0,3))./sqrt(80),squeeze(mean(mean(meanpre(1:6,12,:,freqs>9 & freqs<15),3),4)))
subplot(1,2,2) ; barwitherr(squeeze(std(mean(meanpre(7:12,12,:,freqs>9 & freqs<15),4),0,3))./sqrt(80),squeeze(mean(mean(meanpre(7:12,12,:,freqs>9 & freqs<15),3),4)))




figure, for i=1:25 ; subplot(5,5,i) ; topoplot(squeeze(ica.icawinv(:,i)),ica.chanlocs) ; title(i) ; end




