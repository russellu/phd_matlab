cd c:/shared/badger/russ/topup ; 
ls
fs = dir('res_f_c*') ; 
for i=1:length(fs) ; 
    f = load_untouch_nii(fs(i).name) ;
    resf = reshape(f.img,[1,numel(f.img)]) ; 
    [kc,km] = kmeans(uint8(mat2gray(resf)*255),5) ; 
    kimg = reshape(km,size(f.img)) ; 
    f.img = kimg ; save_untouch_nii(f,['kimg_',num2str(fs(i).name)]) ; 
end














