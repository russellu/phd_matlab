clear all ; close all ; 
cd c:/shared/glaucoma ; ls 

v = VideoReader('20160311_213947.mp4') ;
vframes = read(v) ; 

for i=1:size(vframes,4)
   grays(:,:,i) = rgb2gray(squeeze(vframes(:,:,:,i))) ;  
end

inds1 = 1420:2100 ; 
inds2 = 3350:3920 ; 


