clear all ; close all ; 
cd c:/shared/tests2 ; ls 
subs = dir('*') ; subs(1:2) = [] ; 

for sb = 1%:length(subs) ; 
    %cd(['c:/shared/ute/',subs(sb).name]) ; ls 
    cd c:/shared/tests2
    mask = load_untouch_nii('mederode.nii.gz') ; 
    maskim = mask.img ; 
    outside = bwconncomp(mask.img==0) ; 
    outim = zeros(size(mask.img)) ; 
    inds = outside.PixelIdxList{1} ; 
    outim(inds) = 1 ; 
    intensitylayers = zeros(size(maskim)) ; 
    prevdil = maskim ; 
    for i=1:12 ; 
        dilmaski = imdilate(prevdil,strel(ones(3,3,3))) .* outim ; 
        clayer = (dilmaski - prevdil) > 0 ; 
        intensitylayers(clayer==1) = i ; 
        if i==1
            prevdil = dilmaski ;
        else
            prevdil = (prevdil + dilmaski) > 0 ; 
        end
    end

    mask.img = intensitylayers.*(maskim==0) ; mask.img(:,:,1:10) = 0 ; 
    save_untouch_nii(mask,'intensity_layers.nii.gz') ; 
    skinsurf = (imdilate(mask.img>0,strel(ones(3,3,3))) .* maskim) ;
    mask.img = skinsurf ; save_untouch_nii(mask,'skinsurf.nii.gz') ; 
end