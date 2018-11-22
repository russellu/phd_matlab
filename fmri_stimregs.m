
clear all ; close all ; 
subs = {'sub_alex','sub_charest','sub_esteban','sub_fabio','sub_gab','sub_gabriella','sub_genevieve','sub_gina','sub_jeremie','sub_julie','sub_katrine','sub_lisa'...
        ,'sub_marc','sub_marie','sub_mathieu','sub_maxime','sub_mingham','sub_patricia','sub_po','sub_russell','sub_sunachakan','sub_vincent'} ;
       
    for s=1:length(subs)  ;      
        cd(['c:/shared/allfmris/',subs{s},'/trigs']) ; ls 
        regs = dir('reg_*') ; 
        for r=1:length(regs) ; 
            curr = load_untouch_nii(regs(r).name) ; 
            allregs(s,r,:,:,:) = curr.img ;          
        end
    end
    
stims = [2,3,1,5,6] ;     
cd('c:/shared/regf1/') ; regs = dir('reg_*') ; 
 for r=1:length(regs)  ;
     rf = load_untouch_nii(regs(r).name) ; 
     allrf1(r,:,:,:) = rf.img ; 
  
 end
mf1 = squeeze(mean(allrf1,1)) ;

mregs = squeeze(mean(allregs,1)) ; 
mregs = 1.1-mat2gray((max(max(max(max(mregs))))-mregs));   
 
inds1 = 15:50; inds2=15:60 ; 
for sl=10:25 ; idx = sl; figure
for i=1:length(stims) ; 
    subplot(1,5,i) ;
    %imagesc(squeeze(mean(allregs(:,stims(i),:,:,17),1)),[-20,20]) ;
    plotoverlayIntensity2D(squeeze(mf1(inds1,inds2,idx)),squeeze(mregs(stims(i),inds1,inds2,idx)),squeeze(mregs(stims(i),inds1,inds2,idx))) ;

end   

end

for i=1:66 ; i
    for j=1:66
        for k=1:39
            [p,atab,~] = anova1(squeeze(allregs(:,[1,2],i,j,k)),[],'off') ; 
            sigvol(i,j,k) = atab{2,5} ; 
            
        end
    end
end

    
    