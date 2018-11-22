box = zeros(500,500) ; 
t = 60 ; fr = 85 ; npts = t*fr ; inds = 1:numel(box) ; 
icount = 1 ; 
for i=1:size(box,2):length(inds) - size(box,2)    
    if mod(icount,2) == 0 
        inds(i:i+size(box,2)) = fliplr(inds(i:i+size(box,2))) ; 
    end
    icount = icount + 1 ; 
end

indincr = length(inds)./npts ; 
indpnts = inds(1:indincr:end) ; 
centx = size(box,1)/2 ; centy = size(box,2)/2 ; 
for i=1:length(indpnts)
    [cx,cy] = ind2sub(size(box),indpnts(i)) ; 
    dist = sqrt( (centx-cx).^2 + (centy-cy).^2 ) ; 
    box(cx,cy) = dist ; 
    
end



















