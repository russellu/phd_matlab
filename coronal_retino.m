clear all ; close all  ;
cd C:\shared\glaucoma\restest ; ls 
stims = load('laurenttest_3')  ; stims = stims.stimtimes ; stims = stims(1,:) ; 
nii = load_untouch_nii('Retinotopic_2016_03_17_WIP_sshEPI_1.5mm_TR2650_SENSE_4_1.nii') ; 
nii.img = double(nii.img) ; nimg = nii.img ; 
TR = nii.hdr.dime.pixdim(5) ; 

stimTRs = round(stims./TR) ; 

peakbase = round(2/TR) ; peakstim = round(10/TR) ; 

% get the delay in each voxel:
ts = ([1,2,15,16]) ; 
for i=1:length(ts)
   epochs(i,:,:,:,:) = squeeze(nimg(:,:,:,stimTRs(ts(i))-peakbase:stimTRs(ts(i))+peakstim)) ; 
end

mepochs = squeeze(mean(epochs,1)) ; 
save_nii(make_nii(mepochs),'mepochs.nii.gz') ;




%{
cat = load_untouch_nii('cat2mm.nii.gz') ; 
maxes = zeros(size(mepochs,1),size(mepochs,2),size(mepochs,3)) ; 
for i=1:size(mepochs,1)
    for j=1:size(mepochs,2)
        for k=1:size(mepochs,3)
            maxes(i,j,k) = (find(squeeze(mepochs(i,j,k,6:28))==max(squeeze(mepochs(i,j,k,6:28))),1) + 6)*TR ; 
        end
    end
end
cat.img = maxes ; save_untouch_nii(cat,'maxes.nii.gz') ; 
%}


