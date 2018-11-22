cd c:/shared/asl2/ ; 
subs = dir('nifti*') ;

for s=1%:length(subs) ; 
    cd(['c:/shared/asl2/',subs(s).name]) ;
    ls
    nontags=dir('reg_mc_nontag_asl*') ; 
    tags=dir('reg_mc_tag_asl*') ; 
    for i=1:length(nontags) ; 
        nontagsi = load_untouch_nii(nontags(i).name) ; 
        tagsi = load_untouch_nii(tags(i).name) ; 
        nontagimgs(i,:,:,:,:) = nontagsi.img ; 
        tagimgs(i,:,:,:,:) = tagsi.img ; 
    end
    % get the triggers 
    trignums = [1,3,5] ;
    for i=1:length(trignums)
        trig = load(['stimTimes',num2str(trignums(i)),'.mat']) ;
        trigs{i} = trig.stimTimes ; 
    end
    
    
    asldiffs = tagimgs - nontagimgs ;    
    tlims = [16,40]
    trialbase = tlims(1) ; % seconds of baseline to epoch
    trialtask = tlims(2) ; % seconds of task to epoch
    imgi = squeeze(asldiffs(1,:,:,:,:)) ; 
    boldTR = 8 ;
    nbasevols = round(trialbase/boldTR) ; ntaskvols = round(trialtask/boldTR) ; 
    epochlength = nbasevols + ntaskvols + 1 ;    
    trigsi = trigs{1} ; 
    trigmat = cell2mat(trigsi) ; 
    trigtimesi = trigmat(1:2:end) ; 
    trigtypesi = trigmat(2:2:end) ; 
    onsets = ceil(trigtimesi/boldTR) ; 
    trigtypes = [1,2,6] ; 
    ntrials = 3 ; 
    epochs = zeros(max(size(trigtypes)),ntrials,size(imgi,1),size(imgi,2),size(imgi,3),epochlength) ; 
    stimcounts = ones(1,max(size(trigtypes))) ;
    for onset=1:length(onsets) ;
        currentepoch = imgi(:,:,:,onsets(onset)-nbasevols : onsets(onset)+ntaskvols) ; 
        currenttype = trigtypesi(onset) ;
        currenttypeindex = find(trigtypes==currenttype) ;
        epochs(currenttypeindex,stimcounts(currenttypeindex),:,:,:,:) = currentepoch ; 
        stimcounts(currenttypeindex) = stimcounts(currenttypeindex) + 1 ; 
    end   
    meanepochs = squeeze(mean(mean(epochs,1),2)) ; 
    
    
    
end
