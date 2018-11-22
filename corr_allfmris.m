clear all ; close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 

for sb=1:length(subs)
    cd(['c:/shared/allfmris/sub_',subs{sb}]) ; ls 
    %{
    cd trigs ; 
    stimes = dir('stimTimes*mat') ; 
    for st=1:length(stimes) ; sti = load(stimes(st).name) ; sti = sti.stimTimes ; sti = cell2mat(sti) ; allstims(st,:) = sti(1:2:end) ; end 
    meantimes = mean(allstims,1) ; 
    vec = zeros(1,490) ; for i=1:length(meantimes) ; vec(round(meantimes(i)):round(meantimes(i))+1) = 1 ; end 
    conved = conv(vec,spm_hrf(1),'full') ; conved = conved(1:length(vec)) ; 
    conved = imresize(conved,[1,245]) ; 
    cd .. ; 
    fmri = load_untouch_nii('denoise.nii.gz') ; 
    corrs = voxcorr(fmri.img(:,:,:,20:end-20),conved(20:end-20)) ; 
    f1 = load_untouch_nii('f1.nii.gz') ; 
    f1.img = corrs ; save_untouch_nii(f1,'cleancorrs.nii.gz') ; 
    %}
    cleancorrs = load_untouch_nii('cleancorrs.nii.gz') ; 
    subplot(4,6,sb) ; imagesc(mean(cleancorrs.img,3)) ; 
end
