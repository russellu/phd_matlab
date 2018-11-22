clear all ; close all ; 
cd c:/shared/DICOM_sam/

asl = load_untouch_nii('up.nii.gz') ; 
rest = asl.img(:,:,:,1:80) ; 
stim = asl.img(:,:,:,81:end) ; 

restctrl = stim(:,:,:,1:80) ; 
restlab = stim(:,:,:,81:end) ; 
inds = 1:size(restctrl,4) ; 
for t=1:size(restctrl,4) ; 
    tinds = find(inds >= t-1 & inds <= t+1) ; 
    meant = mean(restlab(:,:,:,tinds),4) ; 
    subs(:,:,:,t) = restctrl(:,:,:,t) - meant ; 

end

ref = load_untouch_nii('stimref.nii.gz') ; 
ref.img = subs ; save_untouch_nii(ref,'stimsub.nii.gz') ; 

stimvec = zeros(1,640) ; 
stimvec(1:154) = 1 ; stimvec(155:155+153) = 0 ; stimvec(154*2:154*2+153) = 1 ; 
hrf = spm_hrf(1) ; 
convstim = conv(stimvec,hrf,'full') ; convstim = convstim(1:length(stimvec)) ; 
resconv = imresize(convstim,[1,80]) ; 
corrs = voxcorr(subs,resconv) ; 
cbfmask = mean(subs,4) > 15 ; 
cbfmask(:,30:end,:) = 0 ;
maskcorrs = corrs.*double(cbfmask) ; 

left = medfilt3(maskcorrs < -.25) ; 
right = medfilt3(maskcorrs > .25) ; 
leftinds = find(left==1) ; 
rightinds = find(right==1) ; 

meanst = load_untouch_nii('meanstim.nii.gz') ; 
meanst.img = left ; save_untouch_nii(meanst,'left_hemi_stim.nii.gz') ; 
meanst.img = right ; save_untouch_nii(meanst,'right_hemi_stim.nii.gz') ; 

rest = load_untouch_nii('restsub.nii.gz') ; 
meanrest = mean(rest.img,4) ; 
meanstim = mean(subs,4) ; 
rstimt = [2:19,40:60] ; 
lstimt = [20:39,61:78] ; 

mright = mean(subs(:,:,:,rstimt),4) ; 
mleft = mean(subs(:,:,:,lstimt),4) ; 

left_restcbf = meanrest(leftinds) ; % resting cbf in left stimulation inds
right_restcbf = meanrest(rightinds) ; % resting cbf in right stimulation inds

left_stimcbf = meanstim(leftinds) ; % both stimulus cbf in left stimulation inds
right_stimcbf = meanstim(rightinds) ; % both stimulus cbf in right stimulation inds

left_leftcbf = mleft(leftinds) ; % left stimulus cbf in left stimulation inds
right_leftcbf = mleft(rightinds) ; % left stimulus cbf in right stimulation inds

left_rightcbf = mright(leftinds) ; % right stimulus cbf in left stimulation inds
right_rightcbf = mright(rightinds) ; % right stimulus cbf in right stimulation inds

clear rest
rest{1} = left_restcbf ; 
rest{2} = right_restcbf ; 
clear leftstim
leftstim{1} = left_leftcbf ;
leftstim{2} = right_leftcbf ; 
clear rightstim
rightstim{1} = left_rightcbf ; 
rightstim{2} = right_rightcbf ; 
save('rest','rest') ; save('leftstim','leftstim') ; save('rightstim','rightstim') ; 


figure, 
subplot(2,2,1) ; 
bar([mean(left_restcbf),mean(right_restcbf)]) ; 
set(gca,'XTickLabel',{'rest, right ctx','rest, left ctx'}) ; title('REST')  ; ylabel('mean CBF (a.u)') ; 
subplot(2,2,2) ; 
bar([mean(left_stimcbf),mean(right_stimcbf)]) ; 
set(gca,'XTickLabel',{'mean hemi, right ctx','mean hemi, left ctx'}) ; title('BOTH stimulation (10min avg)')  ; ylabel('mean CBF (a.u)') ; 
subplot(2,2,3) ; 
bar([mean(left_leftcbf),mean(right_leftcbf)]) ;
set(gca,'XTickLabel',{'left hemi, right ctx','left hemi, left ctx'}) ; title('LEFT stimulation')  ; ylabel('mean CBF (a.u)') ; 
subplot(2,2,4) ; 
bar([mean(left_rightcbf),mean(right_rightcbf)]) ; 
set(gca,'XTickLabel',{'right hemi, right ctx','right hemi, left ctx'}) ; title('RIGHT stimulation') ; ylabel('mean CBF (a.u)') ; 
suptitle('SUM CBF for rest, avg. stim, left hemi stim, right hemi stim') ; 

figure,
subplot(1,2,1) ; 
bar([mean(left_rightcbf)-mean(left_restcbf),mean(right_leftcbf)-mean(right_restcbf)]) ; 
set(gca,'XTickLabel',{'right ctx','left ctx'}) ; title('SUM contra-rest') ; ylabel('mean CBF (a.u)') ; 
subplot(1,2,2) ; 
bar([mean(left_rightcbf)-mean(left_restcbf),mean(right_leftcbf)-mean(right_restcbf)]) ; 
set(gca,'XTickLabel',{'right ctx','left ctx'}) ; title('MEAN contra-rest') ; ylabel('mean CBF (a.u)') ; 
suptitle('DIFFERENCE of SUM and MEAN for right and left cortex') ; 








