cd('C:\shared\badger\EEG-fMRI 2015-10-01\') ; clear all ; close all ; 
stimes = dir('*.mat') ; 
for i=1:length(stimes)
    timesi = load(stimes(i).name) ; 
    timesi = cell2mat(timesi.stimTimes) ; 
    tsecs(i,:) = timesi(1:2:64) ; 
    tstims(i,:) = timesi(2:2:64) ; 
end
vol3 = load_untouch_nii('fMRI_Russell_2015_10_01_WIP_PB-fMRI_SENSE1_5_(saturation)_SENSE_7_1.nii') ; 
TR = 3.244 ; 
hrf = spm_hrf(TR) ;
t3 = tsecs(3,:) ; 
tcourse = zeros(1,size(vol3.img,4)) ; 
stimes = round(t3/TR) ; tcourse(stimes) = 1 ; 
conved = conv(tcourse,hrf,'full') ; conved = conved(1:length(tcourse)) ; 
trimg = vol3.img(:,95:162,:,:) ; 
corrs = voxcorr(trimg,conved) ;
cd trimg ; ls ; a = load('melodic_mix') ; 
comps = load_untouch_nii('melodic_IC.nii.gz') ; comps = comps.img ; 
rgb(:,:,:,1) = uint8(mat2gray(squeeze(comps(:,:,:,30)))*255) ; 
rgb(:,:,:,2) = uint8(mat2gray(squeeze(mean(trimg,4)))*255) ; 
rgb(:,:,:,3) = uint8(mat2gray(squeeze(comps(:,:,:,20)))*255) ; 
for i=1:62 ; figure ; imshow(squeeze(rgb(:,:,i,:))) ; end



%{

a = vol3.img(:,:,:,1) ; 
az = a==0 ; 
trim = a(:,95:162,:) ; 
v3img = vol3.img ; trimg = v3img(:,95:162,:,:) ; 

cvar = squeeze(mean(vol3.img,4))./squeeze(std(double(vol3.img),0,4)) ; 
mn = squeeze(mean(vol3.img,4)) ; 
sdev = squeeze(std(double(vol3.img),0,4)) ; 

rgb = zeros(256,256,62,3) ; 
rgb(:,:,:,1) = uint8(mat2gray(sdev)*255) ; 
rgb(:,:,:,2) = uint8(mat2gray(mn)*255) ; 
rgb(:,:,:,3) = uint8(mat2gray(cvar)*255) ; 
for i=1:62 ; figure ; imshow(squeeze(rgb(:,:,i,:))) ; end
%}
