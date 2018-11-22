function disp3d(img,range)

if nargin == 1 
    maxIM = max(max(max(img))) ;
    minIM = min(min(min(img))) ;
else 
    maxIM = range(2) ; 
    minIM = range(1) ; 
end
    
 

for x=1:100
for i=1:size(img,3)
    imagesc(squeeze(img(:,:,i)),[minIM,maxIM]) ; title(['slice = ',num2str(i)]) ; colormap jet
    getframe ; 
   pause(0.01)
end
end




end