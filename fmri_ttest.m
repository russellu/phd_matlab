clear all ; close all 
subs = {'sub_alex','sub_charest','sub_esteban','sub_fabio','sub_gab','sub_gabriella','sub_genevieve','sub_gina','sub_jeremie','sub_julie','sub_katrine','sub_lisa'...
        ,'sub_marc','sub_marie','sub_mathieu','sub_maxime','sub_mingham','sub_patricia','sub_po','sub_russell','sub_sunachakan','sub_vincent'} ;

    cd C:\shared\ATLASES ; mni = load_nii('MNI152_T1_1mm.nii.gz') ; 
    
for sub = 1:length(subs)
   disp(subs{sub}) ; 
   cd(['c:/shared/allfmris/',subs{sub},'/ants']) ;  
   percs=dir('warp_t1_perc*') ; 
   for perc=1:length(percs) ;
       nii = load_nii(percs(perc).name) ; 
       allpercs(sub,perc,:,:,:) = nii.img ; 
       
       
   end
end
    
% get the t-values in each voxel
 mask = mni.img > mean(mean(mean(mni.img))) ;
ostims = [2,3,4,5,6] ; 
for o=1:length(ostims) ;  
    tbrain = zeros(size(mni.img)) ; 
    for i=1:182 ; disp(i) ;
        for j=1:218
            for k=1:182
                if mask(i,j,k) == 1
                    [h,p,ci,stats] = ttest(squeeze(allpercs(:,1,i,j,k)),squeeze(allpercs(:,ostims(o),i,j,k))) ; 
                    tbrain(i,j,k) = stats.tstat ; 
                end
            end
        end
    end
    allts(o,:,:,:) = tbrain ; 
end
    
   icount = 1 ;  for i=1:2:182 ; subplot(9,12,icount) ; imagesc(rot90(squeeze(allts(1,:,:,i))),[-10,10]) ; icount = icount +1 ;  end

    
    
    