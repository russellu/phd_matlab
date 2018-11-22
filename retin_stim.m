%%%% retinotopic mapping stimulus: 
% minimum 30 second turn-around time for any stimulus location
clear all ; close all ; 
vismask = zeros(1920,1080) ; 
pad = min(size(vismask))/2 ; % two times the padding
icount = 1 ;
while icount < 1000 ; 
    randx = rand ; randy = rand ; 
    rxi = round(randx*(size(vismask,1)-pad)-size(vismask,1)/2 + pad/2) ; 
    ryi = round(randy*(size(vismask,2)-pad)-size(vismask,2)/2 + pad/2) ; 
    fovdist = sqrt(rxi.^2 + ryi.^2) ; 
    rxi_ind = rxi + size(vismask,1)/2 ; 
    ryi_ind = ryi + size(vismask,2)/2 ;  
    circsize = fovdist/4 ; 

    if (sum(sum(vismask(rxi_ind-circsize:rxi_ind+circsize,ryi_ind-circsize:ryi_ind+circsize))) == 0)
        rpts(icount,:) = [rxi,ryi] ; icount = icount + 1 ; 
        [mgx,mgy] = meshgrid(-rxi_ind:size(vismask,1)-rxi_ind-1,-ryi_ind:size(vismask,2)-ryi_ind-1) ; 
        circsize = fovdist/4 ; 
        if circsize<25 ; circsize=25 ; end
        circmask = (sqrt(mgx.^2 + mgy.^2) < (circsize))' ; 
        vismask = vismask + circmask*5 ; 
        vismask(isinf(vismask)) = 0 ; 
        threshmask = (vismask) > mean(mean(vismask)) ; 
        imagesc(threshmask) ; getframe ; pause(.02) ; 
        vismask(vismask>0) = vismask(vismask>0)-.5 ; 
        vismask(vismask<0) = 0 ; 
        vismask = imfilter(vismask,fspecial('gaussian',25,25)) ; 
    end
end
plot(rpts(:,1),rpts(:,2),'o')

