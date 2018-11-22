cd e:/swi_atlas; ls ;

bet_regs = dir('bet_reg*'); 

for betreg=1:length(bet_regs)
swi = load_untouch_nii(bet_regs(betreg).name); 
vein_img = swi.img - imfilter(swi.img,fspecial('gaussian',15,5)); 
[sv,si] = sort(vein_img(:)); 
dilmask = imdilate(swi.img==0, strel(ones(9,9,9))); 

vein_img(si(150000:end)) = 0; 
thresh_veins = double((vein_img<-150) .* (dilmask==0)); 
%disp3d(imfilter(thresh_veins,fspecial('gaussian',15,10)))
swi.img = imfilter(thresh_veins,fspecial('gaussian',15,10));
figure,imagesc(squeeze(swi.img(:,:,50))) ; 

save_untouch_nii(swi,['thresh_',bet_regs(betreg).name]); 

end

