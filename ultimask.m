cd c:/shared/lastute/alex;  ls ; 
rute = load_untouch_nii('gs_20.nii.gz') ; ruteorig = double(rute.img) ; %ruteorig = ruteorig - imfilter(ruteorig,fspecial('gaussian',80,40)); 
z2 = load_untouch_nii('fnirt/mask.nii.gz') ;  
z2.img = medfilt3(imerode(z2.img,strel(ones(3,3,3)))) ; 
layer2 = imdilate(z2.img==2,strel(ones(5,5,5))) ; 
prevdil = double(z2.img==1) ; outim = double(z2.img==0) ; 
intensitylayers = zeros(size(prevdil)) ;
for i=1:12 ; 
    dilmaski = imdilate(prevdil,strel(ones(3,3,3))) .* outim ; 
    clayer = (dilmaski - prevdil) > 0 ; 
    intensitylayers(clayer==1) = i ; 
    if i==1
        prevdil = dilmaski ;
    else
        prevdil = (prevdil + dilmaski) > 0 ; 
    end
    intensitylayers = intensitylayers.*(~layer2) ; 
end
layers = intensitylayers ;










