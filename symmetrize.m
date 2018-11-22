function meanimg = symmetrize(ute)  
% function utesim = symmetrize(ute) ; 
% find the midline of the brain and average mirror images across the
% midline to produce a brain with uniform intensity on either side 
meanimg = zeros(size(ute.img)) ;
fullvis = zeros(size(ute.img)) ; 
for z=1:size(ute.img,3) ; 
    slc = mean(ute.img,3) ; zslice = ute.img(:,:,z) ; 
    slc2 = zeros(size(slc)) ; 
    mthresh = slc > mean(mean(slc)) ; 
    % find the longest continuous slice
    comps = bwconncomp(mthresh) ; compinds = comps.PixelIdxList ; clear lengths ;
    for i=1:length(compinds) ; lengths(i) = length(compinds{i}) ; end ; 
    largest = compinds{find(lengths==max(lengths))} ; 
    largeclust = zeros(size(slc)) ; largeclust(largest) = 1 ; 
    continds = find(sum(largeclust)>5) ; 
    clear cx
    for i=1:length(continds)
        linei = largeclust(:,continds(i)) ; 
        cx(i) = sum((1:length(linei)).*linei')./sum(linei) ;    
    end
    goodinds = find(abs(zscore(cx))<2) ; 
    %plot(continds(goodinds),cx(goodinds),'o') ; 
    fo = fit(continds(goodinds)',cx(goodinds)','poly1') ; 
    % visualize the fit
    cgoods = continds(goodinds) ; 
    visimg = zeros(size(slc)) ; 
    fslope = fo.p1 ; finter = fo.p2 - fslope*min(cgoods) ; 
    for i=1:size(slc,2) 
        midi = round(finter + fslope*i) ; 
        visimg(midi,i) = 1 ; 
    end
    fullvis(:,:,z) = visimg ; 
    for i=1:size(slc,2)
        midlinei = find(visimg(:,i)==1) ; 
        linei = zslice(:,i) ; 
        half1 = linei(1:midlinei) ; half2 = linei(midlinei+1:end) ; 
        linelengths = [length(half1),length(half2)] ; minlength = min(linelengths) ; 
        half2flip = flipud(half2) ; 
        meanhalf = (half1(length(half1)-minlength+1:end)+half2flip(1:minlength)) / 2 ; 
        slc2(midlinei-minlength+1:midlinei,i) = meanhalf ; slc2(midlinei+1:midlinei+minlength,i) = flipud(meanhalf) ; 
    end
    meanimg(:,:,z) = slc2 ; 
end


end









