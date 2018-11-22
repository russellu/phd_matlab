cd('C:\shared\Lyes') ; 

v = VideoReader('stim.mp4') ; 
video = read(v) ;         
for i=1:size(video,4) ;
   grayvid(:,:,i) = rgb2gray(squeeze(video(:,:,:,i))) ; 
    
    
end

%[55:135],[80:220] 

newvid = squeeze(grayvid(55:135,80:220,:)) ; 
boldsecs = 735*0.693 ; 





