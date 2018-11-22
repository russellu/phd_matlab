clear all ; close all 
% script to load and display components from melodic ICA
%cd c:/shared/Test_Russell_2016-03-24/DICOM/melodic ; ls 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
occ = {[65,13],[44,41,28],[60,48],[85,82,67,6],[81,54,36],[77,76,74],[53,49],[66,53,48,38],[77,75,42]} ; 
lingual = {[80,16],[8],[18],[50,5],[60],[53],[31],[6],[16]} ; 
lateral = {[17],[23],[41],[21],[19],[51],[33],[24],[22]} ; 
dmn = {[15,10],[7,2],[12,9,8],[8,1],[1],[33],[43,38],[12,7],[5,3,2]} ; 
motor = {[71,22],[20],[21],[32,11],[53,51,10],[31],[40,5],[88,71,49,8],[51,21,15]} ; 
muscle = {[85],[],[65],[63,61],[85],[87],[40],[60],[81]} ; 
csf = {[5],[5],[3],[31],[48],[6],[29],[33],[29]} ; 
motion = {[33,29],[12],[57],[76,72,30],[59,27],  [34],[63],[81,45],[66,44,43]} ; 
eyes = {[76],[],[70],[80],[26],[55],[55],[51],[41]} ; 
misc_cortex = {[56,54,53,39,35,32,24,9],[34,25,23,15],[30,28,26,31,19,3],[65,55,49,25,22,20,16,14,13,9,4],[45,32,28,20,17,16,15,11,6,5,4],[69,68,60,59,29,25,20,18,17,15,13,12],[33,30,12],[72,61,52,47,42,35,34,80,23,22,14,11],[64,33,32,27,24,13,11,7]} ; 
white = {[51],[10],[27],[19],[8],[36],[20],[17],[30]} ; 
parietal = {[49],[7],[36],[34,12,7],[9],[71,9],[38],[20],[12,6]} ; 
cbllm = {[31],[18],[58],[75],[77],[26],[26],[29],[36]} ; 
back = {[1],[21],[36],[71],[69],[16],[58],[19,15],[18]} ;  


for sub=2 ; 
cd(['c:/shared/newbadger_mri/',subs{sub},'/melodic/']) ; 
all=dir('*') ; all(1:2) = [] ; 
meanNII = load_untouch_nii('mean.nii.gz') ; 
[cx,cy,cz] = centmass3(meanNII.img) ;   
comps = load_untouch_nii('melodic_IC.nii.gz') ;     
axial_anat = squeeze(mean(meanNII.img(:,:,cz-5:cz+5),3)) ; 
sagittal_anat = rot90(squeeze(mean(meanNII.img(cx-5:cx+5,:,:),1))) ; 
mix = load('melodic_mix') ; 
%compinds = {[16,25,14,12,46,64,80],[8,23,2,31,40,41],[18,40,9,15,48,61],[4,21,1,81,84,67],[53,50,14,74,79,80],[17,20,2,41,75,77],[31,33,31,48,50]} ; % fmri

for comp = 1:size(comps.img,4)
    figure ; 
    subplottight(2,2,1) ;
    axial = (mean(squeeze(comps.img(:,:,:,comp)),3)) ; 
    plotoverlayIntensity2D(axial_anat,sqrt(mat2gray(axial)),axial,0) ; 
    subplottight(2,2,2) ; 
    sagittal = rot90(squeeze(mean(squeeze(comps.img(:,:,:,comp)),1))) ; 
    plotoverlayIntensity2D(sagittal_anat,sqrt(mat2gray(sagittal)),sagittal,0) ; 
    subplottight(2,1,2) ;  plot(squeeze(mix(1:735,comp))) ; 
    suptitle(['comp',num2str(comp)]) ; 
end

end