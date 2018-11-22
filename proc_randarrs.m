cd C:\shared\papsaves ; 
ls 

ra = load('randarrs.mat') ; 
ra = ra.randarrs ; 


for i=1:19 ;
    imgi = squeeze(ra(:,:,i)) ; 
    smoothi = imfilter(imgi,fspecial('gaussian',3,1.5)) ; 
    
    
    stds(i) = std(imgi(:)) ; 
    smoothstds(i) = std(smoothi(:)) ; 
    
end

nedges = [1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2] ; 




