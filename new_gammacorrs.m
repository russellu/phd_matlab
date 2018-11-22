subs = {'charest','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','lisa','marc','marie','mathieu','maxime','mingham',...
    'patricia','po','russell','sunachakan','tah','vincent'} ; 

for sub=2%:length(subs) 
    cd(['c:/shared/all_white_normals/a1_good/sub_',subs{sub}]) ;    
    bersp = load('bersp.mat') ; bersp = bersp.bersp ; 
    figure,imagesc(squeeze(mean(mean(mean(bersp([1,3,5],1:2,:,:,:))))),[-3,3]) ; axis xy ; 
    allbers(sub,:,:) = squeeze(mean(mean(mean(bersp([1,3,5],1:2,:,:,:))))) ; 
    
    
    
    
end