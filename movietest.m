clear all ; close all ; 
cd c:/users/butr2901/downloads ; ls 
v = VideoReader('Star Wars- Return of the Jedi VI - Battle of Endor (Space Only) 1080p.mp4') ; 
vframes = read(v) ; 

bwmov = uint8(zeros(size(vframes,1),size(vframes,2),size(vframes,4))) ;
for i=1:size(vframes,4) ; disp(i) ; 
   bwmov(:,:,i) = uint8(rgb2gray(squeeze(vframes(:,:,:,i)))) ; 
   %bwmov(:,:,i) = bwmov(:,:,i) - imfilter(bwmov(:,:,i),fspecial('gaussian',100,5)) ; 
    
end
clear vframes ;



movie = load_untouch_nii('C:\Users\butr2901\Downloads\bwmov.nii') ; 

epoch1 = squeeze(movie.img(:,:,10:10+30*20)) ; 
epoch2 = squeeze(movie.img(:,:,1000:1000+30*20)) ; 
epoch3 = squeeze(movie.img(:,:,2000:2000+30*20)) ; 
epoch4 = squeeze(movie.img(:,:,3000:3000+30*20)) ; 
epoch5 = squeeze(movie.img(:,:,4000:4000+30*20)) ; 
save_nii(make_nii(epoch1),'epoch1.nii') ; 
save_nii(make_nii(epoch2),'epoch2.nii') ; 
save_nii(make_nii(epoch3),'epoch3.nii') ; 
save_nii(make_nii(epoch4),'epoch4.nii') ; 
save_nii(make_nii(epoch5),'epoch5.nii') ; 