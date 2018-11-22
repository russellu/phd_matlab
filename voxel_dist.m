cd c:/shared/ute/MONG_01_RB ; ls 
clear all ; close all ; 
elecorder = {'FP1','FPZ','FP2','AF8','AF4','GND','AF3','AF7','F7','F5','F3','F1','FZ','F2','F4','F6','F8','FT10','FT8','FC6','FC4','FC2','REF','FC1','FC3','FC5','FT7','FT9',...
    'T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP10','TP8','CP6','CP4','CP2','CPZ','CP1','CP3','CP5','TP7','TP9','P7','P5','P3','P1','PZ','P2','P4','P6','P8',...
    'PO8','PO4','POZ','PO3','PO7','O1','OZ','O2'} ; 
coords = load('mricoords.mat') ; coords = coords.mricoords ; 
brain = load_untouch_nii('brain_in_ute.nii.gz') ; 
ebrain = zeros(size(brain.img)) ; 
for i=1:size(coords,2)
    ebrain(coords(1,i),coords(2,i),coords(3,i)) = i ; 
end
brain.img = imdilate(ebrain,strel(ones(5,5,5))) ; save_untouch_nii(brain,'ebrain.nii.gz') ; 

%%% get the new locations from the low resolution electrodes
for i=1:65 ;
    braini = ebrain==i ; 
    [cx,cy,cz] = centmass3(braini==1) ; 
    lowlocs(1,i) = cx ; lowlocs(2,i) = cy ; lowlocs(3,i) = cz ; 
end

lowresbrain = load_untouch_nii('res_brain.nii.gz') ; lowimg = lowresbrain.img ; 
lrb = lowresbrain ; 
binbrain = lowimg > 0 ; 
braininds = find(binbrain==1) ; 
[bx,by,bz] = ind2sub(size(binbrain),braininds) ; 
for i=1:size(coords,2)
    xd(i,:) = lowlocs(1,i)-bx ; 
    yd(i,:) = lowlocs(2,i)-by ; 
    zd(i,:) = lowlocs(3,i)-bz ;  
end
eucls = sqrt(xd.^2 + yd.^2 + zd.^2) ; 
clear xd yd zd ; 
distbrains = zeros(size(binbrain,1),size(binbrain,2),size(binbrain,3),65) ; 
zbrain = zeros(size(binbrain)) ; 
for i=1:size(eucls,1)
    zbrain(braininds) = eucls(i,:) ; 
    distbrains(:,:,:,i) = zbrain ; 
end
clear eucls 
multi = lowresbrain ; multi.img = (mat2gray(distbrains)*10).^2 ; 
multi.hdr.dime.dim = [4,96,96,60,65,1,1,1] ; multi.hdr.dime.pixdim(5) = 1 ; save_untouch_nii(multi,'elecbrain.nii.gz') ; 

distbrains = (mat2gray(distbrains)*10).^2 ; 
wghts = load('wghts') ; wghts = wghts.wghts ; 
rlabs = load('rlabs') ; rlabs = rlabs.rlabs ;  

% correlate the distance in each voxel with the weight map
gnd = find(strcmp('GND',elecorder)) ; 
ref = find(strcmp('REF',elecorder)) ; 
elecorder([gnd,ref]) = [] ; 
distbrains(:,:,:,[gnd,ref]) = [] ;
for i=1:length(elecorder)
   eeginds(i) = find(strcmpi(elecorder{i},rlabs)) ;  
   
end
% eeginds => the indices of the weight file 
corrbrains = zeros(size(distbrains,1),size(distbrains,2),size(distbrains,3),64) ; 
wghtalign = wghts(eeginds,:) ; 
for i=1:64 ; disp(i) ; 
   corrbrains(:,:,:,i) = voxcorr(distbrains,wghtalign(:,i)) ;  
end
corrbrains(isnan(corrbrains)) = 0 ; 
multi = lowresbrain ; multi.img = corrbrains ; 
multi.hdr.dime.dim = [4,96,96,60,64,1,1,1] ; multi.hdr.dime.pixdim(5) = 1 ; save_untouch_nii(multi,'corrbrain.nii.gz') ; 


