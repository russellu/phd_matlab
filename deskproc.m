cd c:/
desk = imread('desk_bright (1).jpg') ;
graydesk = rgb2gray(desk) ; 
dthresh = graydesk < 220 ;
invthresh = graydesk > 220 ;
cd c:/newdesk
for g=50:5:220
for i=1:3
    newdesk(:,:,i) = medfilt2(double(squeeze(desk(:,:,i))).*double(dthresh) + invthresh.*g,[3,3]);
    
end
imwrite(uint8(newdesk),['graylvl_',num2str(g),'_','desk1.png']) ;
end
%imshow(uint8(newdesk))
