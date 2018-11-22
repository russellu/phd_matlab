cd c:/shared/Douleur_2016-04-28/ ; ls 

nozoom = load_untouch_nii('bp_desp_topup_nozoom.nii.gz') ; 
zoom = load_untouch_nii('bp_desp_topup_zoom.nii.gz') ; 

tr = 1750 ; 
secs = 1:500*1.75 ; 
ts = zeros(1,round(500*1.75)) ; 


for i=1:6
   ts((i-1)*120+20*(i-1)+20 : (i-1)*120+20*(i-1)+20+120) =1 ;  
end
hrf = spm_hrf(1.75) ; 
tstr = imresize(ts,[1,500]) ; 
conved = conv(tstr,hrf,'full') ; conved = conved(1:500) ; 

corrzoom = voxcorr(zoom.img,conved) ; 
corrnozoom = voxcorr(nozoom.img,conved) ; 
avgcorrs = (corrzoom + corrnozoom) / 2 ; 

f = load_untouch_nii('f_zoom.nii.gz') ; 
f.img = avgcorrs ; save_untouch_nii(f,'corrs.nii.gz') ; 