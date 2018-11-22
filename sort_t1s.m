subs = {'alex','biz','cloud','dave','dina','felix','genevieve','jeremie','karl','nic','pierre','russell','sukh','tegan','terry','valerie'} ; 
for sb=1%:length(subs)
    cd(['c:/shared/simul_fs/',subs{sb},'/mri']) ; 
    t1 = load_untouch_nii('T1.nii.gz') ; 
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,2,1) ; imagesc(squeeze(t1.img(128,:,:))) ; 
    subplot(2,2,2) ; imagesc(squeeze(t1.img(:,128,:))) ; 
    subplot(2,2,3) ; imagesc(squeeze(t1.img(:,:,128))) ; 
end
