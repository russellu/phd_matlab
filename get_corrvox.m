clear all ; close all ; 
subs= {'alex','dina','genevieve','jeremie','karl','russell','tegan','valerie'} ;

for sub=1:length(subs) 
sub = subs{sub} ; 
cd(['c:/shared/badger_mri/',sub,'/nii']) ; 

hrf1 = load('gamma_hrf_1.txt') ; hrf2 = load('gamma_hrf_2.txt') ; 

bps=dir('bp*gamma*') ; 

bp1 = load_untouch_nii(bps(1).name) ; 
bp2 = load_untouch_nii(bps(2).name) ; 

corrs1 = voxcorr(bp1.img(:,:,:,50:end-50),hrf1(50:end-50)) ; 
corrs2 = voxcorr(bp2.img(:,:,:,50:end-50),hrf2(50:end-50)) ; 

mcorrs = (corrs1 + corrs2)./2 ; 

g = find(mcorrs>.1) ; 
[sv,si] = sort(mcorrs(g),'descend') ; 
zcorrs = zeros(size(mcorrs)) ; 
zcorrs(g(si(1:500))) = 1 ; 
zcorrs = imdilate(zcorrs,strel(ones(3,3,3))) ; 
bw = bwconncomp(zcorrs) ; 
pix = bw.PixelIdxList ; 
sizes = cellfun(@length,pix) ; 

maxsz = find(sizes==max(sizes)) ; 
maxinds  = pix{maxsz} ; 
zclust = zeros(size(mcorrs)) ; 
zclust(maxinds) = 1 ; 

fref = load_untouch_nii('fref.nii.gz') ; 
fref.img = mcorrs ; 
save_untouch_nii(fref,'meancorrs.nii.gz') ; 
fref.img = zclust ; 
save_untouch_nii(fref,'zclust.nii.gz') ; 

end

