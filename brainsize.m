clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 
%subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'} ; 
for sb=1:length(subs); disp(sb) ; 
    cd(['c:/shared/allfmris/sub_',subs{sb}]) ;  
   % cd(['c:/shared/badger_struct/',subs{sb}]) ;  

    brain = load_untouch_nii('gm.nii.gz') ; 
    sbrain = brain.img>0 ; 
    sumbrains(sb) = sum(sbrain(:)) ;
    
    
end
