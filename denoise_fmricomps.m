% automatic ICA denoising of FMRI time series (removing bad components and
% saving their indices as a text file so that melodic can subtract them

clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie','tegan'} ; 

compinds = {[16,25,14,12,46,64,80],[8,23,2,31,40,41],[18,40,9,15,48,61],[4,21,1,81,84,67],[37,18,2,32,41,54,79],[53,50,14,74,79,80],[17,20,2,41,75,77],[31,33,31,48,50]} ; 

for s=1:length(subs)
cd (['c:/shared/badger_mri/',subs{s},'/nii/melodic']) ; ls 
mix = load('melodic_mix') ; 
meanf = load_untouch_nii('mean.nii.gz') ; 
[cx,cy,cz] = centmass3(meanf.img) ; 
anat = squeeze(mean(meanf.img(:,:,cz-6:cz+6),3)) ; 
comps = load_untouch_nii('melodic_IC.nii.gz') ; 
%{
figure
for i=1:length(compinds{s}) ; 
    subplottight(1,3,i) ;
    clear img
    img(:,:) = mat2gray(squeeze(max(comps.img(:,:,cz-6:cz+6,compinds{s}(i)),[],3))) ;
    %img(:,:,2) = mat2gray(squeeze(mean(comps.img(:,:,14:20,i),3))) ;
    %img(:,:,3) = mat2gray(squeeze(mean(comps.img(:,:,21:28,i),3))) ;
    plotoverlayIntensity2D(anat,mat2gray(log(abs(img))),img,0) ; title(i) ; 
    set(gca,'XTickLabel',[],'YTickLabel',[]) ; 
    
end
%}



figure
for i=1:size(comps.img,4) ; 
    subplottight(9,10,i) ;
    clear img
    img(:,:) = mat2gray(squeeze(max(comps.img(:,:,cz-6:cz+6,i),[],3))) ;
    %img(:,:,2) = mat2gray(squeeze(mean(comps.img(:,:,14:20,i),3))) ;
    %img(:,:,3) = mat2gray(squeeze(mean(comps.img(:,:,21:28,i),3))) ;
    plotoverlayIntensity2D(anat,mat2gray(log(abs(img))),img,0) ; title(i) ; 
    set(gca,'XTickLabel',[],'YTickLabel',[]) ; 
    
end

end

