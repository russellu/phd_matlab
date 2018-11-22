clear all ; close all ; 
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
    ls
    s = load('posepochs') ; s = s.posepochs ;
    sizevox(sub) = size(s,2) ; 
    bcorrs = (s(:,:,:,:) - repmat(mean(s(:,:,:,1:2),4),[1,1,1,11]))./repmat(mean(s(:,:,:,1:2),4),[1,1,1,11]) ; 
    alls(sub,:,:,:,:) = squeeze(mean(bcorrs,2)) ; 
    
end




