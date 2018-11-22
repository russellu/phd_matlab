cd c:/shared/gammaregs ; ls 

clocs = dir('colorlocs*gz') ; 
count =1 ; 
for c=1:length(clocs) ; 
   cloci = load_untouch_nii(clocs(c).name) ; 
   inds = unique(cloci.img(find(cloci.img ~=0))) ;
   if length(inds)==65
      for i=1:65 ; 
          [cx(count,i),cy(count,i),cz(count,i)] = centmass3(cloci.img==i) ; 
      end
      count = count +1 ; 
   end
end

meancoords = round([mean(cx);mean(cy);mean(cz)]) ; 
meanelecbrain = zeros(size(cloci.img)) ; 
for i=1:size(meancoords,2)
   meanelecbrain(meancoords(1,i),meancoords(2,i),meancoords(3,i)) = i ;  
end
meanelecbrain = imdilate(meanelecbrain,strel(ones(5,5,5))) ; 
cloci.img = meanelecbrain ; 
save_untouch_nii(cloci,'meanelecbrain.nii.gz') ; 
cloci.img = meanelecbrain > 0; 
save_untouch_nii(cloci,'binmeanelecbrain.nii.gz') ; 
