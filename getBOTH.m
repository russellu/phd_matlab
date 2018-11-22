
function [eeg,fmri,icawinv,comps,rawcomps,subjects] = getBOTH() 

subjects = {
    'alex'
    'charest'
    'esteban'
    'fabio'
    'gab'
    'gabriella'
    'genevieve'
    'gina'
  %  'guillaume'
    'jeremie'
    'julie'
    'katrine'
    'lisa'
    'marc'
    'marie'
    'mathieu'
    'maxime'
    'mingham'
    'patricia'
    'po'
    'russell'
    'sunachakan'
  %  'tah'
  %  'thititip'
    'vincent'
} ; 


for sub=1:size(subjects,1) ; 
    cd(['c:/shared/allfmris/sub_',subjects{sub}]) ;
    s = load('posepochs') ; s = s.posepochs ;
    sizevox(sub) = size(s,2) ; 
    bcorrs = (s(:,:,:,:) - repmat(mean(s(:,:,:,1:2),4),[1,1,1,11]))./repmat(mean(s(:,:,:,1:2),4),[1,1,1,11]) ; 
    fmri(sub,:,:,:,:) = squeeze(mean(bcorrs,2)) ; 
    
end

fmri  = squeeze(mean(fmri,3)) ; 

for sub=1:size(subjects,1) ; 
    disp(subjects{sub}) ; 
    cd(['c:/shared/allres/',subjects{sub}]) ;
    s = load('compeeg') ; 
    s = s.compeeg ; 
    eeg(sub,:,:,:) = s ;   
    winv = load('icawinv') ; icawinv(sub,:,:) = winv.icawinv ; 
    inds = load('inds') ; comps(sub,:) = inds.inds ; 
    raw = load('rawcomps') ; rawcomps(sub,:,:,:,:) = raw.rawcomps ; 
end

end



