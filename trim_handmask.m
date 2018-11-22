clear all ; close all ; 
cd c:/shared/lastute/ ;
mongs=dir('*') ; mongs(1:2) = [] ; 
mricoordn = 1 ; 
for m=1%:length(mongs) ; 
    disp(mongs(m).name) ; 
    cd(['C:\shared\lastute\',mongs(m).name,'\ants\']) ;
    mask = load_untouch_nii('utemask.nii.gz') ; 
    resute = load_untouch_nii('res_ute.nii.gz') ; 
    gs = resute.img - imfilter(resute.img,fspecial('gaussian',60,30)) ; 

    shell = mask.img - imerode(mask.img,strel(ones(3,3,3))) ; 
    
    shellinds = find(shell~=0 ) ; 
    vals = resute.img(shellinds) ; 
    [kc,km] = kmeanscustom(uint8(mat2gray(vals)*255),2) ; 
    zs = zeros(size(shell)) ; 
    zs(shellinds) = km ; 
end