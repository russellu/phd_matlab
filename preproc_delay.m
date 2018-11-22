subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie'} ; 

for s=1:length(subs) 
    cd(['c:/shared/badger_mri/',subs{s},'/nii/']) ; ls   
    nii = load_untouch_nii('delay_mean.nii.gz') ; 
    img = squeeze(nii.img(:,:,:,3)) ;   
    [~,km] = kmeans(uint8(mat2gray(reshape(img,[1,numel(img)]))*255),3) ;  
    kimg = reshape(km,size(img)) ; 
    figure,
    for i=1:33 ; subplot(3,11,i)
        imagesc(squeeze(kimg(:,:,i)),[0,4]) ;
    end
    anat = load_untouch_nii('f_topup_mc_retino_allstims_01.nii.gz') ; 
    dilimg = kimg==3 ; 
    dilimg = imdilate(dilimg,strel(ones(3,3,3))) ; 
    bw = bwconncomp(dilimg) ; 
    plist = bw.PixelIdxList ; 
    lengths = cellfun('length',plist) ; 
    clustimg = zeros(size(dilimg)) ; 
    clustimg(plist{lengths==max(lengths)}) = 1 ; 
    clustimg = medfilt3(clustimg) ; 
    anat.img = clustimg ; 
    
    save_untouch_nii(anat,'kmask.nii.gz') ; 
    
    
end
