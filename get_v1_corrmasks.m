%clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 


for sub=1:length(subs)
%%%% BOLD FMRI processing:
cd(['c:/shared/badger_mri/',subs{sub},'/nii/']) ; 
mask = load_untouch_nii('meancorrs.nii.gz') ; 
medmask = medfilt3(mask.img>.1) ; 
bw = bwconncomp(medmask) ; 
pix = bw.PixelIdxList ; 
lens = cellfun(@length,pix) ; 
maxlen = find(lens==max(lens)) ; 
clust1 = pix{maxlen} ; 
zs = zeros(size(medmask)) ; 
zs(clust1) = 1 ; 
mask.img = zs ; save_untouch_nii(mask,'epi_v1_cluster.nii.gz') ; 
allsums(sub) = sum(sum(sum(mask.img))) ; 
end

clear corrs ; 
for i=1:50  ;
    for j=1:3
       corrs(i,j) = corr2(squeeze(mbersp(:,j,i)),allsums') ;   
    end
end

for sub=1:length(subs)
    cd c:/shared/v1_and_modulation
    mkdir(subs{sub}) ; 
    cd(subs{sub}) ; 
    subbersp = squeeze(mbersp(sub,:,:)) ; 
    save('baseline_corrected_eeg','subbersp') ; 
end