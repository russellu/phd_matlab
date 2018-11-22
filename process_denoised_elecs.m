cd c:/shared/lastute/alex ; 
ls 
for i=1:8
   nii = load_untouch_nii(['colorimg',num2str(i),'.nii.gz']) ;  
   img = nii.img ; 
   coordsi = load(['mricoords_',num2str(i)]) ; 
   allcoords(i,:,:) = coordsi.mricoords ; 
   for j=1:65
       eleci = img==j ; 
       [cx,cy,cz] = centmass3(eleci) ; 
       coords(i,j,:) = [cx,cy,cz] ;
   end
end

plot(squeeze(mean(std(coords,0,1),3))) ; hold on ;
plot(squeeze(mean(std(allcoords,0,1),2)),'r') ; 






