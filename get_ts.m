clear all ; close all ; 
cd c:/shared/simdenoise/russell/ ; ls 
eeg = pop_loadbv('.','retino_gamma_01.vhdr') ; 

labs = {eeg.urevent.type} ; 
lats = {eeg.urevent.latency} ; lats = cell2mat(lats) ; 
r128s = find(strcmp('R128',labs)) ; 
stims1 = find(strcmp({'S  1'},labs)) ; stims2 = find(strcmp({'S  2'},labs)) ; stims3 = find(strcmp({'S  3'},labs)) ; 
allstims = [stims1,stims2,stims3] ; 
stims = sort(lats(allstims))-lats(r128s(1)) ; 
stimtrs = round((stims/5000)/ 0.693) ; 
zvec = zeros(1,735) ; 

for i=1:length(stimtrs) ; zvec(stimtrs(i):stimtrs(i)+round(5/0.693)) = 1 ; end
conved = conv(zvec,spm_hrf(0.693),'full') ; conved = conved(1:length(zvec)) ; 
dlmwrite('stim_hrf_1.txt',conved') ; 

%nii = load_untouch_nii('c:/shared/ignatius_mri/russell/nii/bp_reg_topup_mc_retino_gamma_01.nii.gz') ; 
%corrs = voxcorr(nii.img,conved) ; 
%dlmwrite('stim_hrf_1.txt',conved) ; 







