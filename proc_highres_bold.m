cd('c:/shared/3D fMRI 2017-01-12/DICOM') ; 

lyes = load_untouch_nii('IM_0042.nii.gz') ; 
img = lyes.img ; 
ts = zeros(1,120) ; ts(61:end) = 1 ; 
ts = repmat(ts,[1,5]) ; 
newtrs = round(600*(1/1.15)) ; 
ts = imresize(ts,[1,newtrs]) ; 
ts = conv(ts,spm_hrf(1.15),'full') ; 
ts = ts(1:size(img,4)) ; 
corrs = voxcorr(img,ts) ; 
ref = load_untouch_nii('lyescat.nii.gz') ; 
ref.img = corrs ; save_untouch_nii(ref,'lyes_corrs.nii.gz') ; 


cd('c:/shared/3D fMRI 2017-01-12/DICOM') ; 
bold = load_untouch_nii('IM_0033.nii.gz') ; 
mbold = bold.img(:,:,:,67:end) + bold.img(:,:,:,34:66); 
ts = zeros(1,120) ; ts(61:end) = 1 ; 
ts = repmat(ts,[1,5]) ; 
newtrs = ceil(600*(1/19))+1 ; 
ts = imresize(ts,[1,newtrs]) ; 
corrs = voxcorr(mbold,ts) ; 
ref = load_untouch_nii('ref3d.nii.gz') ; 
ref.img = corrs ; save_untouch_nii(ref,'corrs3d.nii.gz') ; 









