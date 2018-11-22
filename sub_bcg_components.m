clear all ; close all; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'}; 
scans = {'bp_reg_topup_mc_retino_allstims_01','bp_reg_topup_mc_retino_allstims_02','bp_reg_topup_mc_retino_gamma_01','bp_reg_topup_mc_retino_gamma_02','bp_reg_topup_mc_retino_movie','bp_reg_topup_mc_retino_rest'}; 
eegscans = {'retino*allstim*01*set','retino*allstim*02*set','retino*gamma*01*set','retino*gamma*02*set','retino*movie*set','retino*rest*set'}; 

minpeaks = [180,160,170,200,200,180,180,170,140]; % min inter-peak distance (samples) for BCG findpeaks (adjust per subject)
comps = [1,4,3,3,2,3,2,3,3]; % BCG ICA component (set after visualizing results in line 24)
polarities = [1,1,-1,-1,-1,-1,-1,-1,1]; % polarity of component (peaks positive (1) or negative (-1)?) (set after visualizing)
epochls = {[100,200-15],[100-20,200-35],[100,200-35],[100-25,200],[100,200],[100,200],[100,200-10],[100,200],[100-10,200-60]}; % samples to take around each bcg peak 

for sb=1:length(subs)
cd(['E:\badger_eeg\',subs{sb}]);

dats = dir('retino*set'); % brainvision output,
for dat=1:length(dats) % for all saved datasets

EEG = pop_loadset(dats(dat).name); % load using header
% replace BCG artifact with noise (so it doesn't impact ICA), this channel 
% will differ between experimental setups
EEG.data(32,:) = rand(1,length(EEG.data))*5;

bdat = eegfiltfft(EEG.data,EEG.srate,1,128); % for average subtraction (we are not removing values under 1Hz)
hdat = eegfiltfft(EEG.data,EEG.srate,3,128); % for individual BCG epoch correlation (more accurate to exclude low-freq)

allweights = load('allweights.mat'); allweights = allweights.allweights; % weights from get_bcg_components
acts = squeeze(allweights(:,:))*bdat; 
template = acts(comps(sb),:)*polarities(sb) ; % get the bcg component time series selected in previous script
[pks,locs] = findpeaks(smooth(template),'MINPEAKDISTANCE',minpeaks(sb)); % get the peaks, use subject specific value for minpeaks


elength = epochls{sb}(1) + epochls{sb}(2); % length of one bcg epoch (samples)
eprev = epochls{sb}(1); epost = epochls{sb}(2); % before and after peak number of samples
eps = zeros(64,length(locs),elength+1); % lower highpass epochs
heps = zeros(64,length(locs),elength+1); % higher highpass epochs
epinds = zeros(length(locs),elength+1); % epoch indices
for i=2:length(locs)-2 % get the BCG epochs using locs from findpeaks
    eps(:,i,:) = bdat(:,locs(i)-eprev:locs(i)+epost); 
    heps(:,i,:) = hdat(:,locs(i)-eprev:locs(i)+epost); 
    epinds(i,:) = locs(i)-eprev:locs(i)+epost;
end
meaneps = squeeze(mean(eps,2)); % mean lower highpass epochs (all channels, mean across epoch)


[bsv,bsi] = sort(sum(abs(meaneps),2),'descend'); % sort channels by epoch magnitude 
clear corrs
for i=1:size(heps,1) % create correlation matrix by correlating all epochs
    corrs(i,:,:) = corr(squeeze(heps(i,:,:))');  

end
mcorrs = squeeze(mean(corrs(:,:,:),1)); % mean across epochs, can use bsi to only average channels with more BCG
[sv,si] = sort(mcorrs,2,'descend') ; % sort according to highest correlation for each epoch
neweeg = EEG; % create new EEG for subtracted data
for i=2:size(epinds,1)- 2
    for j=1:64
        neweeg.data(j,epinds(i,:)) = squeeze(EEG.data(j,epinds(i,:))) - squeeze(mean(eps(j,si(i,5:40),:),2))' ; % subtract mean of top correlating (5-40) epochs
    end
end
% visualize some stuff (template BCG component, all channels' epochs, topopot of BCG epoch strength in each channel)
%figure,subplot(2,2,1); plot(template); vline(locs);
%subplot(2,2,2); plot(meaneps'); 
%sums = sum(abs(diff(neweeg.data,1,2)),2); 
%allsums(dat,:) = sums; 
%subplot(2,2,3); topoplot(sums,neweeg.chanlocs,'electrodes','numbers')
%alleps(dat,:,:) = meaneps; 
% save subtracted set, merge for further analysis
pop_saveset(neweeg,['bcg_',strrep(dats(dat).name,'.dat','')]); 
if dat==1; merged = neweeg; else merged = pop_mergeset(merged,neweeg); end
if dat==1; rawmerged = EEG; else rawmerged = pop_mergeset(rawmerged,EEG); end

end

figure,
%subplot(1,2,1) ; plot(squeeze(mean(alleps,1))') ; title(subs{sb}); 
plot(template); vline(locs); title(subs{sb}); 
end





% run ICA on merged data and visualize components
%{
figure,topoplot(mean(allsums,1),neweeg.chanlocs,'electrodes','numbers');
filtnew = eegfiltfft(merged.data,EEG.srate,3,128); 
[weights,sphere] = runica(filtnew(:,1:3:end),'maxsteps',128) ; 
winv = pinv(weights*sphere); 
figure,for i=1:64; subplottight(5,13,i) ; topoplot(winv(:,i),EEG.chanlocs) ; end
acts = weights*sphere*neweeg.data; 
%}


