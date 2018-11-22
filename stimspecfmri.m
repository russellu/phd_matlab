%stimspecfmri.m requires that epoch_fmri is run previous.
% gets the stimulus specific activations maps. 
% needs the variable fimgs to be set


clear allstims ; 
    tcounts = ones(6,1); 

for trial=1:size(allparams,2) ;
    for s=1:size(trtimes,1)
        allstims(allparams(1,trial,s),tcounts(allparams(1,trial,s)),:,:,:,:) = fimgs(trial,:,:,:,trtimes(s)-2:trtimes(s)+4) ; 
        tcounts(allparams(1,trial,s)) = tcounts(allparams(1,trial,s)) + 1 ; 


    end
end




tvals = (squeeze(mean(mean(allstims(:,:,:,:,:,5:7),2),6)-mean(mean(allstims(:,:,:,:,:,1:2),2),6)))./squeeze(std(mean(allstims(:,:,:,:,:,1:2),6),0,2)./sqrt(45)) ;
save('tvals','tvals') ; 

