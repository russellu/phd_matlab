clear all; close all
cd e:/fmris/badger_alex; 


ls 
f1 = load_untouch_nii('f1.nii.gz'); 
mask = f1.img>mean(f1.img(:)); 
resmask = mask(:); maskinds = find(resmask==1); 
nii = load_untouch_nii('reg_topup_mc_retino_gamma_01.nii.gz'); 
resimg = reshape(nii.img,[numel(nii.img(:,:,:,1)),size(nii.img,4)]); 

[weights,sphere] = runica(resimg(maskinds,:)','maxsteps',128); 

winv = pinv(weights*sphere); 

acts = weights*sphere*resimg(maskinds,:)'; 


newacts = zeros(size(resimg)); 



cd e:/badger_eeg/


nc = 100; 
zica = fastICA(resimg(maskinds,:),nc,'negentropy');
zicimg = zeros(nc,numel(nii.img(:,:,:,1))); 
for i=1:size(zica,1)
    zicimg(i,maskinds) = zica(i,:); 
end
reszic = reshape(zicimg',[size(nii.img(:,:,:,1)),nc]); 
for i=1:100 ; subplottight(10,10,i) ; imagesc(rot90(squeeze(max(reszic(:,:,:,i),[],1)))) ;colormap jet; end


for i=1:size(acts,1)
    newacts(maskinds,i) = acts(i,:); 
end

resacts = reshape(newacts,size(nii.img)); 
nii.img = resacts; save_untouch_nii(nii,'resacts.nii.gz'); 














