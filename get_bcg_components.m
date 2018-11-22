clear all ; close all; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'}; 
scans = {'bp_reg_topup_mc_retino_allstims_01','bp_reg_topup_mc_retino_allstims_02','bp_reg_topup_mc_retino_gamma_01','bp_reg_topup_mc_retino_gamma_02','bp_reg_topup_mc_retino_movie','bp_reg_topup_mc_retino_rest'}; 
eegscans = {'retino*allstim*01*set','retino*allstim*02*set','retino*gamma*01*set','retino*gamma*02*set','retino*movie*set','retino*rest*set'}; 

comps = [1,4,3,3,2,3,2,3,3]; % BCG ICA component (set after visualizing results in line 24)
polarities = [1,1,-1,-1,-1,-1,-1,-1,1]; % polarity of component (peaks positive (1) or negative (-1)?) (set after visualizing)

for sb=1%:length(subs)
cd(['E:\badger_eeg\',subs{sb}]);
dats = dir('*set'); % brainvision output,
for dat=1:length(dats)%  get all and merge
EEG = pop_loadset(dats(dat).name); %load using header
EEG.data(32,:) = rand(1,length(EEG.data))*5; 
if dat==1; merged = EEG; else merged = pop_mergeset(EEG,merged); end 
end

filt = eegfiltfft(merged.data,merged.srate,1,128); % high pass
[weights,sphere] = runica(filt(:,1:4:end),'maxsteps',128) ; % ICA
allweights = weights*sphere; % save weights
acts = weights*sphere*filt; 
save('allweights','allweights'); 

% visualize top 12 components
figure,for i=1:12 ; subplot(3,4,i); plot(squeeze(acts(i,10000:40000))); title(i); end; suptitle(subs{sb}); 
end


