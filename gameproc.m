cd c:/game/alexAttackingProc/ ; ls 
pngs = dir('*png') ; 
for png=1:length(pngs)
   imgs(:,:,:,png) = imread(pngs(png).name) ;  
end

whites = imgs ; whites(imgs==0) = 255 ; 
for i=1:14 ; subplot(4,4,i) ; imshow(squeeze(whites(:,:,:,i))) ; end


