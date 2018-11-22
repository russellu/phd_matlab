clear all ; close all ; 

subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 
freqs = 1:1:80 ; 
badrestts = {[1:35,420:450],...
             [1:35,420:450],...
             [1:35,35:40,338:348,420:450],...
             [1:35,135:155,245:255,420:450],...
             [1:35,18:60,170:180,270:280,375:385,420:450],...
             [1:35,48:70,110:120,160:170,193:203,220:236,299:306,332:341,414:421,420:450],...
             [1:35,144:160,420:450]} ; 
         
badmoviets = {[1:35,41:49,700:734],...
              [1:35,109:124,200:221,444:455,700:734],...
              [1:35,170:210,464:485,600:615,630:640,677:690,700:734],...
              [1:35,120:140,310:316,480:495,540:560,584:590,695:708,700:734],...
              [1:35,350:360,495:510,570:580,700:734],...
              [1:35,21:38,74:80,110:118,156:166,176:185,220:245,295:305,340:355,380:386,410:418,429:435,458:465,495:505,514:525,585:600,662:675,700:734],...
              [1:35,275:290,590:600,650:662,670:690,700:734]} ;


sub = 1 ; 
cd(['c:/shared/badger_eeg2/',subs{sub}]) 
bcgs = dir('bcgica*movie*set') ;

for b=1%:length(bcgs)
    EEG = pop_loadset(bcgs(b).name) ; 
    EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
    acts = EEG.icaweights*EEG.icasphere*EEG.data ; 
    facts = zeros(length(freqs),size(acts,1),size(acts,2)) ; 
    winv = pinv(EEG.icaweights*EEG.icasphere) ; 
    %figure,for i=1:32 ; subplot(4,8,i) ; topoplot(squeeze(winv(:,i)),EEG.chanlocs) ; title(i) ; end
    clear facts 
    for i=1:length(freqs)
        facts(i,:,:) = eegfiltfft(acts,250,freqs(i)-1,freqs(i)+1) ; 
    end
    [s,f] = spectopo(acts,0,250,'plot','off') ; 
    
    events = {EEG.urevent.type} ; 
    times = cell2mat({EEG.urevent.latency}) ; 
    trigs = find(strcmp('R128',events)) ; 
    tlats = times(trigs) ; 
    pow = facts(:,:,tlats(1):tlats(end)).^2 ; 
    for i=1:size(pow,1)
        for j=1:size(pow,2)
            pow(i,j,:) = imfilter(squeeze(pow(i,j,:)),fspecial('gaussian',[100,1],25)) ; 
        end
    end

    convFactor = 250*0.693 ; % number of samples per TR     
    resSize = size(pow,3)/convFactor ; 
    ntrs = round(resSize) ; 
    respow = zeros(size(pow,1),size(pow,2),resSize) ; 
    for i=1:size(pow,1) ; disp(i) ; 
        for j=1:size(pow,2)
            respow(i,j,:) = imfilter(imresize(squeeze(pow(i,j,:)),[ntrs,1]),fspecial('gaussian',[12,1],2)) ; 
        end
    end
    cd(['c:/shared/newbadger_mri/',subs{sub}]) ; ls
    bp = load_untouch_nii('bp_reg_topup_mc_retino_movie.nii.gz') ; bpimg = bp.img ; 

    corr = load_untouch_nii('corrs.nii.gz') ; 
    mask = medfilt3(corr.img > .25) ; inds = find(mask==1) ; [i1,i2,i3] = ind2sub(size(mask),inds) ; 
    for i=1:length(inds)
        ts(i,:) = squeeze(bpimg(i1(i),i2(i),i3(i),:)) ; 
    end
    mts = mean(ts,1) ; mts = mts(1:ntrs) ; 
    clear zrestts 
    zpower = zeros(1,size(respow,3)) ; zrestts(badmoviets{sub}) = 1 ; resclusts = bwconncomp(zrestts) ;   

    newrespow = respow ;  
    for i=1:length(resclusts.PixelIdxList)
        clusti = resclusts.PixelIdxList{i} ; 
        for j=1:size(respow,1)
            for k=1:size(respow,2)
                minclust = clusti(1) ; maxclust = clusti(end) ; 
                cluststep =  ((respow(j,k,maxclust)-respow(j,k,minclust)) / length(clusti)) ; 
                if cluststep ~= 0
                    linevals = respow(j,k,minclust):cluststep:respow(j,k,maxclust) ; 
                    linevals = linevals(1:length(clusti)) ; 
                    newrespow(j,k,clusti) = linevals ; 
                else
                    linevals = ones(1,length(clusti)) * minclust ; 
                    newrespow(j,k,clusti) = linevals ; 
                end            
            end           
        end   
    end
    respow = newrespow ; 
    inds = 50:10:size(respow,3)-50 ; 
    clear cmat ; 
    for idx=1:length(inds)
        for i=1:size(respow,1)
            for j=1:64
                tsi = mts(inds(idx)-35:inds(idx)+35) ; 
                respi = respow(i,j,inds(idx)-35:inds(idx)+35) ; 
                cmat(idx,i,j,:) = xcorr(tsi,respi,25,'coeff') ; 
            end
        end
    end
    mcmat = squeeze(mean(cmat,1)) ;
    chani = squeeze(cmat(:,:,3,:)) ;
end




