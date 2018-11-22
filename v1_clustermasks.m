clear all ; close all ; 
subs= {'alex','dina','genevieve','jeremie','karl','russell','tegan','valerie'} ;

for s=1:length(subs) ; 
    clear ts tsinds corrinds maxcorrs
    cd(['c:/shared/badger_mri/',subs{s},'/nii']) ; ls 
    retinds = load_untouch_nii('maxinds_t1.nii.gz') ;                        
    retcorrs = load_untouch_nii('meancorrs_t1.nii.gz') ; 
    
    corrmask = retcorrs.img > .15 ; 
    
    maxinds = find(retcorrs.img>.15) ; 
    [sv,si] = sort(retcorrs.img(maxinds),'descend') ; 
    seedclust = maxinds(si(1:2000)) ; 
    
    corrimg = zeros(size(retcorrs.img)) ; 
    corrimg(seedclust) = 1 ;  
    
    
    dilmask = imdilate(corrimg,strel(ones(21,21,21))).*corrmask ; 
    
    %imagesc(max(retcorrs.img,[],3)) ; 
    medmask = medfilt3(retcorrs.img > .15) ; 
    maskedinds = medmask.*retinds.img ; 
    binmasks = zeros(size(medmask)) ; 
    icount = 1 ; 
    for i=1:10:80
        binmasks(maskedinds>=i & maskedinds<i+10) = icount  ;
        icount = icount + 1 ; 
    end
    binmasks = binmasks.*dilmask ; %((imdilate(binmasks,strel(ones(3,3,3))).* double(binmasks==0)) + binmasks ) ; 
    
    for i=1:max(binmasks(:))
        retinds.img = medfilt3(binmasks==i) ; disp(i) ; 
        save_untouch_nii(retinds,['mask_',num2str(i),'_clustinds.nii.gz']) ;   
    end
    
    
    
    retinds.img = binmasks ; save_untouch_nii(retinds,'clustinds.nii.gz') ; 
    
    %{
    [cx,cy,cz] = ind2sub(size(binmasks),find(binmasks>0)) ; 
    anglevals = binmasks(binmasks>0) ; 
    roivec = [cx,cy,cz,anglevals] ; 
    [kc,km] = kmeans(roivec,20) ; 
    
    binmasks(find(binmasks>0)) = kc ; 
        retinds.img = binmasks ; save_untouch_nii(retinds,'clustinds.nii.gz') ; 
    %}
    
    
end




