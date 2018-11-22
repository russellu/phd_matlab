cd c:/shared/badger/felix/ica_reg_mc_fmri_1.nii.gz
maps = load_untouch_nii('melodic_IC.nii.gz') ; 
meanf = load_untouch_nii('mean.nii.gz') ; 
mapims = maps.img ; 
mixes = load('melodic_mix') ; 
badzs = find(max(zscore(mixes,1)) > 6 | min(zscore(mixes,1)) < -6) ; 
goodzs = find(max(zscore(mixes,1)) < 6 & min(zscore(mixes,1)) > -6) ; 
gmix = mixes(:,goodzs) ; 

for i=1:size(gmix,2) ; figure,plot(gmix(:,i)) ; end






