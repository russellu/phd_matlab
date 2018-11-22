%cd c:/shared/SWI_TOF_DIFF/DICOM; ls 
clear all ; close all ;
cd('C:\shared\fMRI+stimulation sujet 2 2015-07-16') ; ls 
%mri = load_untouch_nii('im_0010.nii.gz')  ;
mri = load_untouch_nii('Test_Francis_2015-07-16_WIP_fMRI_stim2_SENSE_15_1.nii')  ;
mri = double(mri.img) ; 
brain = zeros(size(mri,1),size(mri,2),size(mri,3),(floor(size(mri,4)./2))) ; 
for i=1:size(mri,1) ; disp(i) ; 
    for j=1:size(mri,2);
        for k=1:size(mri,3)
            f = abs(fft(squeeze(mri(i,j,k,:)))) ; 
            brain(i,j,k,:) = f(1:floor(size(mri,4)./2)) ;       
        end
    end
end
save_nii(make_nii(brain),'brain.nii.gz') ; 