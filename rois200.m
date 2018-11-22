clear all ; close all ; 
rois = load_untouch_nii('c:/shared/ATLASES/kbrain.nii.gz') ; 
cd c:/shared/claudie ; ls 
all=dir('*') ; all(1:2) = [] ; 
for file = 1:length(all) 
    cd(['c:/shared/claudie/',all(file).name]) ; 
    prenii = load_untouch_nii('bp_mni_denoised_pre.nii.gz') ; pre_img = prenii.img ; 
    postnii = load_untouch_nii('bp_mni_denoised_post.nii.gz') ; post_img = postnii.img ; 
    
    clear pre_roi post_roi
    for i=1:200
        inds = find(rois.img==i) ;
        [x,y,z] = ind2sub(size(rois.img),inds) ; 
        for j=1:length(x)
           pre_roi(j,:) = squeeze(pre_img(x(j),y(j),z(j),:)) ;  
           post_roi(j,:) = squeeze(post_img(x(j),y(j),z(j),:)) ;      
        end
        pre_all(i,:) = mean(pre_roi) ; 
        post_all(i,:) = mean(post_roi) ; 
    end
    
    for i=1:200 ; disp(i) ; 
        corrpre(i,:,:,:) = voxcorr(pre_img,pre_all(i,:)) ; 
        corrpost(i,:,:,:) = voxcorr(post_img,post_all(i,:)) ; 
        
    end
    
    cat = load_untouch_nii('mnicat.nii.gz') ; 
    cat.img = permute(corrpre,[2,3,4,1]) ; 
    save_untouch_nii(cat,'corrpre.nii.gz') ; 
    cat.img = permute(corrpost,[2,3,4,1]) ; 
    save_untouch_nii(cat,'corrpost.nii.gz') ;     
    
    
end




