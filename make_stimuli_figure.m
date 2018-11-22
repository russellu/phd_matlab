clear all ; close all; 

% make the retinotopic stimuli:
[xg,yg] = meshgrid(-400:400,-400:400) ; 
width = length(xg) ; height = length(yg) ; 
circmask = double(sqrt(xg.^2 + yg.^2) < 350);
fullfield = sqrt(xg.^2 + yg.^2) < 150 ; fullfield = fullfield.*circmask ;
left = (xg<0) .* circmask ; right = (xg>=0) .* circmask ; top = (yg<=0).*circmask ; bot = (yg>0).*circmask ; 
[th,rh] = cart2pol(xg,yg) ; 
wedge1 = (th>0 & th < pi/4).*circmask ; wedge2 = (th >= pi/4 & th < pi/2).*circmask ; 
wedge3 = (th>=pi/2 & th < 3*pi/4).*circmask ; wedge4 = (th>=3*pi/4 & th <= pi).*circmask ; 
wedge5 = (th<=0 & th < -(pi-pi/4)).*circmask ; wedge6 = (th >= -(pi-pi/4) & th < -(pi-pi/2)).*circmask ;  
wedge7 = (th >= -(pi-pi/2) & th < -(pi-3*pi/4)).*circmask ; wedge8 = (th >= -(pi-3*pi/4) & th <= 0).*circmask ; 
allwedge(:,:,1) = wedge1 ; allwedge(:,:,2) = wedge2 ; allwedge(:,:,3) = wedge3 ; allwedge(:,:,4) = wedge4 ; 
allwedge(:,:,5) = wedge5 ; allwedge(:,:,6) = wedge6 ; allwedge(:,:,7) = wedge7 ; allwedge(:,:,8) = wedge8 ; 
ring1 = double(sqrt(xg.^2 + yg.^2) < 400) & double(sqrt(xg.^2 + yg.^2) >=200) ;
ring2 = double(sqrt(xg.^2 + yg.^2) < 200) & double(sqrt(xg.^2 + yg.^2) >=50) ;
ring3 = double(sqrt(xg.^2 + yg.^2) < 50) & double(sqrt(xg.^2 + yg.^2) >=0) ;
allrings(:,:,1) = ring1 ; allrings(:,:,2) = ring2 ; allrings(:,:,3) = ring3 ; 
quad1 = (th>0 & th<pi/2).*circmask ; quad2 = (th>=pi/2 & th<=pi).*circmask ; 
quad3 = (th>=-pi & th<-pi/2).*circmask ; quad4 = (th>=-pi/2 & th<=0).*circmask ; 
allquads(:,:,1) = quad1  ; allquads(:,:,2) = quad2 ; allquads(:,:,3) = quad3 ; allquads(:,:,4) = quad4 ; 
rots = 22.5:22.5:361 ; 
bar1 = (yg<10 & yg>-10).*circmask ; bar2 = imrotate(bar1,rots(1),'nearest','crop').*circmask ; 
bar3 = imrotate(bar1,rots(2),'nearest','crop').*circmask ; bar4 = imrotate(bar1,rots(3),'nearest','crop').*circmask ;
bar5 = imrotate(bar1,rots(4),'nearest','crop').*circmask ; bar6 = imrotate(bar1,rots(5),'nearest','crop').*circmask ;
bar7 = imrotate(bar1,rots(6),'nearest','crop').*circmask ; bar8 = imrotate(bar1,rots(7),'nearest','crop').*circmask ; 
allbars(:,:,1) = bar1 ; allbars(:,:,2) = bar2 ; allbars(:,:,3) = bar3 ; allbars(:,:,4) = bar4 ; 
allbars(:,:,5) = bar5 ; allbars(:,:,6) = bar6 ; allbars(:,:,7) = bar7 ; allbars(:,:,8) = bar8 ; 
clear allstims ; 
allstims(:,:,1) = fullfield ; allstims(:,:,2) = left ; allstims(:,:,3) = right ; allstims(:,:,4) = top ; allstims(:,:,5) = bot ; 
allstims(:,:,6:13) = allwedge ; allstims(:,:,14:17) = allquads ; allstims(:,:,18:20) = allrings ; 
for i=1:8
   allstims(:,:,i+20) = imrotate(allwedge(:,:,i),22.5,'crop') ;  
end
circmask = double(sqrt(xg.^2 + yg.^2) < 70);
clear rhosin

[xg,yg] = meshgrid(-400:400,-400:400) ; 
[th,rh] = cart2pol(xg,yg) ;
rhosin = sin(rh/4) ;  
lines = sin(xg/4); 


titles = {'full','left','right','top','bottom','bottom right','bottom left','top left','top right','periphery','fovea'};
retstims = [1,2,3,4,5,14,15,16,17,19,20]; 
figure,
for i=1:length(retstims)  
    subplot(1,11,i) ; 
    maski = imresize(squeeze(allstims(:,:,retstims(i))).*allstims(:,:,1),[801,801]) ; 
    allmasks(i,:,:) = maski; 
    fimgs(:,:,i) = (uint8(mat2gray(rhosin.*maski+1)*255)) ;
    limgs(:,:,i) = (uint8(mat2gray(lines.*maski+1)*255)) ;
    finimg = (uint8(mat2gray(rhosin.*maski+1)*255)); 
    ret_stims(:,:,i) = finimg(250:end-250,250:end-250); 
    imshow(finimg(250:end-250,250:end-250)); 
    colormap gray; set(gca,'XTickLabel',[],'YTickLabel',[]);  
end
cd e:\saved
stims_discr = ret_stims ; save('stims_discr','stims_discr'); 

allmasks = allmasks > 0.5; 
figure,
rotangles = [45,90,135,180,225,270,315,360]; 
rotlabels = {'0\circ','45\circ','90\circ','135\circ','180\circ','225\circ','270\circ','315\circ'};
for i=1:length(rotangles)
wedge = imrotate(squeeze(fimgs(:,:,6)),rotangles(i),'bicubic','crop');
subplottight(8,1,i); 
imshow(wedge(250:end-250,250:end-250)); 
stims_wedge(:,:,i) = wedge(250:end-250,250:end-250); 
text(320,100,rotlabels{i}); 
end
save('stims_wedge','stims_wedge'); 

figure,
rotangles = [90,135,180,225,270,315,360,-45]; 
rotlabels = {'0\circ','45\circ','90\circ','135\circ','180\circ','225\circ','270\circ','315\circ'};
for i=1:length(rotangles)
bars = imrotate(squeeze(limgs(:,:,1)),rotangles(i),'bicubic','crop');
subplottight(8,1,i); 
imshow(bars(250:end-250,250:end-250)); 
stims_bars(:,:,i) = bars(250:end-250,250:end-250); 
text(320,100,rotlabels{i}); 
end
save('stims_bars','stims_bars'); 

contrasts = [1,0.3,0.15,0.05,0.01];

for i=1:length(contrasts)
    subplot(1,5,i); 
maski = imresize(squeeze(allstims(:,:,retstims(1))).*allstims(:,:,1),[801,801]) ; 

contimg = (mat2gray(limgs(:,:,1))-.5)*2;
contimg = uint8(contimg*contrasts(i)*127.5+127.5); 
contimg = double(contimg).*double(maski); contimg(maski==0) = 127.5; 
imshow(uint8(contimg(250:end-250,250:end-250))); 
%limgs(:,:,i) = (uint8(mat2gray(lines.*maski+1)*255)) ;
%imshow(limgs(:,:,i)); 
colormap gray; set(gca,'XTickLabel',[],'YTickLabel',[]);  
stims_contrasts(:,:,i) = uint8(contimg(250:end-250,250:end-250));

end
save('stims_contrasts','stims_contrasts'); 






