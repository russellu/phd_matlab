cd C:\Users\butr2901\Desktop\tuning\surfaces

names = {'full','left','right','top','bot','bot_right','bot_left','top_right','top_left','periphery','fovea'}; 

for n=1:length(names)
   imgs(n,:,:,:) = imread(['discrete_',names{n},'.png']); 
    
end

resimgs = reshape(imgs,[size(imgs,1)*size(imgs,2)*size(imgs,3),3]); 
sumres = sum(resimgs,2); for i=1:3; resimgs(sumres==0,i) = 255; end; 
whiteimgs = reshape(resimgs,size(imgs)); 
imagesc(squeeze(whiteimgs(1,200:600,200:800,:)));

clipimgs = whiteimgs(:,200:600,200:800,:); 
cd E:\saved
save('clipimgs','clipimgs'); 