clear all ; close all ; 
cd c:/shared/glaucoma/laurent ; ls 
both = load('laurent_stim_both_delay.mat') ; 
both = both.stimtimes ; 
both = both(1,:) ; 

cd func/DICOM ; 
nii = load_untouch_nii('mc_timing_unwarped.nii.gz') ; 
ntrs = size(nii.img,4) ; 
TR = 1.24 ; 
timing = load_untouch_nii('f1_timing.nii.gz')  ;


% make a time series in ms of the same length as the BOLD
ms = 1:ntrs*TR*1000 ; 
bothms = both*1000 ; 
zms = zeros(1,length(ms)) ; 
for i=1:length(bothms)
   vals = abs(ms-bothms(i)) ; 
   minval = find(vals==min(vals)) ; 
   zms(minval:minval+5000) = 1 ; 
end

stim = imresize(zms,[1,ntrs]) ; 

hrf = spm_hrf(1.24) ; 
conved = conv(stim,hrf,'full') ; conved = conved(1:length(stim)) ; 
corrbrain = voxcorr(nii.img,conved) ; 
timing.img = corrbrain ; 
save_untouch_nii(timing,'corrbrain.nii.gz') ; 
dlmwrite('stim.txt',stim') ; 




