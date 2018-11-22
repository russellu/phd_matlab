rbrain = mbrain; 
midgamma = double((squeeze(mean(rbrain(:,:,:,60:80),4))) < -.15) ; 
%midgamma = (squeeze(mean(rbrain(:,:,:,40:110),4))<-2) ; 
beta = double((squeeze(mean(rbrain(:,:,:,15:25),4))) < -.2) ; 
highgamma = double(squeeze(mean(rbrain(:,:,:,80:110),4))>.15) ; 
clear img xim;
img(:,:,:,1) = highgamma ; img(:,:,:,2) = midgamma ;img(:,:,:,3) = beta ; %img = mat2gray(img) ;  
for i=1:3 ; img(:,:,:,i) = medfilt3(squeeze(img(:,:,:,i)),[3,3,3]) ; end
incr = 3 ; icount = 1;
for i=10:22
    subplot(3,4,icount) ; 
    for j=1:3
        xim(:,:,j) = rot90(squeeze(mean(img(:,:,i:i+incr,j),3))) ;
    end
    
    imshow(rot90(squeeze(mask(:,:,i))), 'InitialMag', 'fit') ;  hold on ; 
    %h = imshow(xim) ;
    im(:,:,1) = rot90((squeeze(macs(:,:,i))>0).*(squeeze(macs(:,:,i))))*2 ; 
    im(:,:,3) = abs(rot90(squeeze(macs(:,:,i)<0).*(squeeze(macs(:,:,i)))))*2 ; 
    h = imshow(im) ; 
    %I = mat2gray(squeeze(mean(xim>0,3))>0) ;
   % I = rot90(abs(squeeze(macs(:,:,i)))*2) ; 
    set(h, 'AlphaData', I)
    title(i) ; 
    icount = icount + 1 ; 
end