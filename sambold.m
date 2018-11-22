cd c:/shared/samtest_asl_bold
bold = load_untouch_nii('bp_mc_lightson.nii.gz') ; 
secs = 350*0.693 ; 
ts = zeros(1,ceil(secs)) ; 
ts(31:60) = 1  ; ts(91:120) = 1 ; ts(151:180) = 1 ; 
conved = conv(ts,spm_hrf(1),'full') ; conved = conved(1:length(ts)) ; 
conved = imresize(conved,[1,size(bold.img,4)]) ; 
corrs = voxcorr(bold.img,conved) ; 

for i=1:33 ; 
    subplottight(3,11,i) ; 
    imagesc(rot90(squeeze(corrs(:,:,i))),[-.75,.75]) ; 
    set(gca,'XTick',[],'YTick',[]) ; 
end













