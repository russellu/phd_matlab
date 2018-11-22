cd C:\shared\Lyes\lyes_nifti ; ls 

% get the stimulus
cd lyes_retinotopic ; ls 
fullfield = load('lyes_fullfield.mat') ; fullfieldtimes = fullfield.stimtimes ; 
retinotopy = load('lyesretinotopy.mat') ; retinotopytimes = retinotopy.stimtimes ; cd .. ; 


% sequence = downright,downleft,upright,upleft,right,left,down,up,rotcw,rotccw,inwards,outwards

% get the full field convolution for the delay image 
fullfield = load_untouch_nii('reg_bp_full_field.nii.gz') ;
ntrs = size(fullfield.img,4) ; 
TR = 1.24 ; 
TRtimes = round(fullfieldtimes./TR) ; 
ideal = zeros(1,ntrs) ; 
for i=1:length(TRtimes) ; ideal(TRtimes(i):TRtimes(i)+round(5/TR)) = 1 ; end
dlmwrite('ideal.txt',ideal') ; 
convideal = conv(ideal,spm_hrf(TR),'full') ; convideal = convideal(1:ntrs) ; 
corrs = voxcorr(fullfield.img(:,:,:,40:end-40),convideal(40:end-40)) ; 
ref = load_untouch_nii('f1_fullfield.nii.gz') ; ref.img = corrs ; save_untouch_nii(ref,'corrs.nii.gz') ; 

















