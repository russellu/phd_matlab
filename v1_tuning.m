%clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie'} ; 

for s=1:length(subs) 
    cd(['c:/shared/badger_mri/',subs{s},'/nii/']) ; ls   
    %nii = load_untouch_nii('delay_mean.nii.gz') ; 
    %img = squeeze(nii.img(:,:,:,3)) ;   
    gamma1 = load_untouch_nii('bp_reg_topup_mc_retino_gamma_01.nii.gz') ; 
    gamma2 = load_untouch_nii('bp_reg_topup_mc_retino_gamma_02.nii.gz') ; 
    conv1= dlmread('gamma_hrf_1.txt') ; 
    conv2= dlmread('gamma_hrf_2.txt') ; 
    
    corrs1 = voxcorr(gamma1.img,conv1) ; 
    corrs2 = voxcorr(gamma2.img,conv2) ; 
    
    mcorrs = (corrs1+corrs2)/2 ; 
    %figure,
    %for i=1:33 ; subplot(4,9,i) ; imagesc(squeeze(mcorrs(:,:,i))>.3,[-1,1]) ; end
    angles = load_untouch_nii('allangles.nii.gz') ; 
    angleimg = angles.img ; 
    
    stdangles = angleimg./repmat(std(angleimg,0,4),[1,1,1,size(angleimg,4)]) ; 
    [ix,iy,iz] = ind2sub(size(mcorrs),find(mcorrs>.3)) ; 
    clear v1
    for i=1:length(ix)
        v1(i,:) = squeeze(stdangles(ix(i),iy(i),iz(i),:)) ; 
    end
    
    [~,si] = sort(mean(v1(:,1:90),2),'descend') ; 
    
    normv1 = v1.*(v1>0) ; 
    normv1(isnan(normv1)) = 0 ; 
    figure,bar(mean(normv1,1))
    allmeans(s,:) = mean(normv1,1) ; 
    
    
 
end
