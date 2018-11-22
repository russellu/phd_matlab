[x,y] = meshgrid(-300:299,-300:299) ; 
[theta,rho] = cart2pol(x,y) ; 
%%%% you want it to expand more quickly as it gets further and further out.
% also ,you want it to get larger and larger as it gets further and further
% out. 
% first thing you need is to get the radii of each annulus (inner and
% outer), just have a function for this.
% oh and also you need to have enough so that it makes a smooth transition
% in time, ie, if the framerate is 60fps you need 600 frames

nframes = 200 ; incr = log(200)/nframes ; 
outer = exp(1:incr:(log(200)+1)) ; % good
% for the inner...same thing, just with a small offset and padded with
% 0s...
fovthresh = 25 ; 
rhos = zeros(600,600,200) ; 
for i=1:length(outer)
    if i < fovthresh+1
        rhos(:,:,i) = rho < outer(i) ; 
    else
        rhos(:,:,i) = rho < outer(i) & rho > outer(i-fovthresh) ;     
    end
end
disp3d(rhos)

