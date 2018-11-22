clear all ; close all ; 
subs = {'sub_lyes','sub_amal','sub_lyndis','sub_valerie','sub_reihaneh','sub_samuel','sub_cesar','sub_angelina','sub_russell','sub_esteban','sub_raphael','sub_greg','sub_felix','sub_pascal'};
for sb=1:length(subs) 
    cd(['e:\orientation_retinotopy\',subs{sb}])
    
    ls 
    
    elecimg= load_untouch_nii('epi_elecs.nii.gz');
    
    for i=1:65
      [cx(sb,i),cy(sb,i),cz(sb,i)] = centmass3(elecimg.img==i);   
    end
    
    retimg = load_untouch_nii('native_rets.nii.gz'); 
    
    retimg.img(:,1:40,:,:) = 0; 
    subplot(4,5,sb) ; imagesc(squeeze(mean(mean(retimg.img,3),4))); 
    
    for i=1:11
       [mx(sb,i),my(sb,i),mz(sb,i)] = centmass3(retimg.img(:,:,:,i)); 
        
    end
 
end
cd E:\nimg_pool\saved
elecorder = load('elecorder.mat') ; elecorder = elecorder.elecorder; 
for i=1:11
    for j=1:14
        for k=58:65
           dists(i,j,k) = norm([mx(j,i) - cx(j,k), my(j,i) - cy(j,k), mz(j,i) - cz(j,k)]);
        end
        
    end
end
mdists = squeeze(mean(dists(:,:,58:65),3)); 


barwitherr(squeeze(std(mdists,0,2))/sqrt(14),mean(mdists,2));



