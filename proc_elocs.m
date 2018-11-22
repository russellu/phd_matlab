clear all ; close all ; 
EEG = pop_loadbv('C:\shared\MONG_01_RB\','MONG_01_RB_FIX_BOX.vhdr') ;
EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
rlabs = {EEG.chanlocs.labels} ; 
cd c:/shared/raw ; allelecs = dir('*MONG*') ; 
for el=1:length(allelecs) ; disp(allelecs(el).name) ; 
cd(['c:/shared/raw/',allelecs(el).name]) ;
cmass = load('cmass.mat') ; cmass = cmass.cmass ; 
cx = cmass{1} ; cy = cmass{2} ; cz = cmass{3} ; 
wm = load_untouch_nii('t1_in_ute.nii.gz') ; wimg = wm.img ; 
rimg = zeros(size(wimg)) ; 
ss = load_untouch_nii('ss_ute.nii.gz') ; ssimg = ss.img>0 ; 
alldists = zeros(size(rimg,1),size(rimg,2),size(rimg,3),63) ; 
for i=1:64 ; disp(i) ; 
[xn,yn,zn] = ndgrid(-cx(i):size(wimg,1)-cx(i)-1,-cy(i):size(wimg,2)-cy(i)-1,-cz(i):size(wimg,3)-cz(i)-1) ; 
distxyz = sqrt(xn.^2+yn.^2+zn.^2) ; 
alldists(:,:,:,i) = distxyz.*double(ssimg) ; 
%rimg = rimg + (sqrt(xn.^2 + yn.^2 + zn.^2) < 35) ; %rimg = max(max(max(rimg)))-rimg ; 
end
save_nii(make_nii(uint8(alldists)),'alldists.nii.gz');
%wm.img = imfilter(rimg,fspecial('gaussian',35,35))  ; save_untouch_nii(wm,'elecsum.nii.gz') ; 
%mask = mat2gray(squeeze(rimg(:,:,cz(i))))/2 ;
%figure,plotoverlayIntensity2D(squeeze(wimg(:,:,cz(i))),mask,squeeze(rimg(:,:,cz(i))),270) ; 
%wm.img = sqrt(rimg) ; save_untouch_nii(wm,'rimg.nii.gz') ;
end


b = load_untouch_nii('sumbrain.nii.gz') ; 
ut1 = load_untouch_nii('t1_in_ute.nii.gz') ; 

for i=1:5:size(ut1.img,3) ; figure,
plotoverlayIntensity2D(squeeze(ut1.img(:,:,i)),mat2gray(squeeze(b.img(:,:,i))),squeeze(b.img(:,:,i)),0) ; 
end




%{
cd c:/shared/delrecs/ ; skull = load_untouch_nii('skfill.nii') ; 
bw2 = imfill(skull.img>0,'holes') ;
er = bw2 ; 

for i=1:8
    diler = imdilate(er,strel(ones(3,3,3))) ; 
    layers(i,:,:,:) = (diler-er)*i ; 
    er = diler ; 
    
end
skull.img = squeeze(sum(layers,1)) ; save_untouch_nii(skull,'layers.nii.gz') ; 
%}
