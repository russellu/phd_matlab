function pupildata = interpolate_blinks(pupildata,threshold,smoothfac)

% interpolate blinks
zmask = pupildata<threshold; zmask = imdilate(zmask,strel(ones(1,smoothfac))); 
bw = bwconncomp(zmask); 
for i=1:bw.NumObjects
   inds = bw.PixelIdxList{i}; 
   if i==1
       val = mean(pupildata([inds(end)+1]));
       pupildata(inds) = val; 
   elseif i==(bw.NumObjects)
       val = mean(pupildata([inds(1)-1]));
       pupildata(inds) = val; 
   else
       val = mean(pupildata([inds(end)+1,inds(1)-1]));
       pupildata(inds) = val; 
   end
end

end