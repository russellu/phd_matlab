function disp4d(input,n)
%%% display a 4d on n subplots
%%% assumes time to be the fourth dimension and z the 3rd

plotStep = floor(size(input,3)/n)  
plotsqr = ceil(sqrt(n)) ; 

for reps = 1:100

for t = 1:size(input,4)

    plotCount = 1 ;
    maxt = max(max(max(squeeze(input(:,:,:,t))))) ; 
    mint = min(min(min(squeeze(input(:,:,:,t))))) ; 
    
    for i=1:plotStep:size(input,3)-plotStep
        subplot(plotsqr,plotsqr,plotCount) ;
        plotCount = plotCount + 1 ; 
        imagesc(squeeze(input(:,:,i,t))) ;
    end
    getframe ; 
end

end

end