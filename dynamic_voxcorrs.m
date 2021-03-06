clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 
compinds = {[16,25,14,12,46,64,80],[8,23,2,31,40,41],[18,40,9,15,48,61],[4,21,1,81,84,67],[53,50,14,74,79,80],[17,20,2,41,75,77],[31,33,31,48,50]} ; % fmri
badcomps = {[1,2,3,6,16,20,29,35:44,46:49,51:64],[1,2,4,5,7,10,11,16,28:64],[1,2,3,4,7,8,14,15,16,17,19,20,22,24,29,30,33,34:39,41:64],... % eeg
    [1,2,3,5,6,10:13,17:19,21:23,25:64],[1:6,9,10,14,16,20:25,27,28,30:33,35:38,41:64],[1:7,9,12,13,15,19:22,24,25,27:31,33,34,36:64],...
    [1:3,5,7,8,9,15,17,21,23,24:30,32,33,35:64]} ;

goodcomps = {[7,8,10,12,14,15,18,19,22,23,28,30],[3,8,12,13,18,25,26],[6,9,11,12,13,18,23,28,31,40],[7,8,20,24],[13,17,18,26,29,39],[8,17,18,23,32,35],[6,11,14,18,20,31,34]} ;
elecorder = {'FP1','FPZ','FP2','AF8','AF4','GND','AF3','AF7','F7','F5','F3','F1','FZ','F2','F4','F6','F8','FT10','FT8','FC6','FC4','FC2','REF','FC1','FC3','FC5','FT7','FT9',...
    'T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP10','TP8','CP6','CP4','CP2','CPZ','CP1','CP3','CP5','TP7','TP9','P7','P5','P3','P1','PZ','P2','P4','P6','P8',...
    'PO8','PO4','POZ','PO3','PO7','O1','OZ','O2'} ; 
elecs = [59,45,31,46,60,9,20,18] ; 

badrestts = {[1:35,420:450],...
             [1:35,420:450],...
             [1:35,35:40,338:348,420:450],...
             [1:35,135:155,245:255,420:450],...
             [1:35,18:60,170:180,270:280,375:385,420:450],...
             [1:35,48:70,110:120,160:170,193:203,220:236,299:306,332:341,414:421,420:450],...
             [1:35,144:160,420:450]} ; 
         
badmoviets = {[1:35,41:49,700:735],...
              [1:35,109:124,200:221,444:455,700:735],...
              [1:35,170:210,464:485,600:615,630:640,677:690,700:735],...
              [1:35,120:140,310:316,480:495,540:560,584:590,695:708,700:735],...
              [1:35,350:360,495:510,570:580,700:735],...
              [1:35,21:38,74:80,110:118,156:166,176:185,220:245,295:305,340:355,380:386,410:418,429:435,458:465,495:505,514:525,585:600,662:675,700:735],...
              [1:35,275:290,590:600,650:662,670:690,700:735]} ;

for sub=1:length(subs)
    
    %%%% BOLD FMRI processing:
    cd(['c:/shared/badger_mri/',subs{sub},'/nii/melodic']) ; 
    mix = load('melodic_mix') ; 
    weights = load_untouch_nii('melodic_IC.nii.gz') ; 
    allweights{sub} = weights.img ; 
    meansub = load_untouch_nii('mean.nii.gz') ; 
    allmeans{sub} = meansub.img ; 
    mixinds = 1:735:735*5 ; 
    for i=1:length(mixinds) ; boldinds{i} = mixinds(i):mixinds(i)+734 ; end
    boldinds{6} = 735*5+1:735*5+1+449 ; 
    newboldinds{1} = boldinds{5} ; newboldinds{2} = boldinds{6} ; 
    boldinds = newboldinds ; 
    clear segmix
    for i=1:length(boldinds) ; smix = eegfiltfft(mix(boldinds{i},:)',1/0.693,0.01,2) ; segmix{i} = smix' ; end
    allsegmix{sub} = segmix ; 
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    restpow = load('post_restpow') ; restpow = restpow.restpow ; 
    moviepow = load('post_moviepow') ; moviepow = moviepow.moviepow ; 
    
    % find the motion from the auto-motion corrected points. (each TR
    % binary)
    restset = dir('*rest*Pulse*set') ; 
    movieset = dir('*movie*Pulse*set') ; 
    
    restmotion = pop_loadset(restset(1).name) ; 
    moviemotion = pop_loadset(movieset(1).name) ; 
   
    figure, 
    subplot(1,2,1) ; 
    bar(squeeze(mean(mean(restpow(elecs,35:end,:),1),2))) ; vline(badrestts{sub}) ; 
    subplot(1,2,2) ;
    bar(squeeze(mean(mean(moviepow(elecs,35:end,:),1),2))) ; vline(badmoviets{sub}) ; title(subs{sub}) ; 
    
    zmoviets = zeros(1,length(moviepow)) ; zmoviets(badmoviets{sub}) = 1 ; movieclusts = bwconncomp(zmoviets) ; 
    zrestts = zeros(1,length(restpow)) ; zrestts(badrestts{sub}) = 1 ; restclusts = bwconncomp(zrestts) ;   
    
    newrest = restpow ;  
    for i=1:length(restclusts.PixelIdxList)
        clusti = restclusts.PixelIdxList{i} ; 
        % interpolate the power at each frequency
        for j=1:size(restpow,2)
            for k=1:size(restpow,1)
                minclust = clusti(1) ; maxclust = clusti(end) ; 
                cluststep =  ((restpow(k,j,maxclust)-restpow(k,j,minclust)) / length(clusti)) ; 
                if cluststep ~= 0
                    linevals = restpow(k,j,minclust):cluststep:restpow(k,j,maxclust) ; 
                    linevals = linevals(1:length(clusti)) ; 
                    newrest(k,j,clusti) = linevals ; 
                else
                    linevals = ones(1,length(clusti)) * minclust ; 
                    newrest(k,j,clusti) = linevals ; 
                end            
            end           
        end   
    end
    
    newmovie = moviepow ;  
    for i=1:length(movieclusts.PixelIdxList)
        clusti = movieclusts.PixelIdxList{i} ; 
        % interpolate the power at each frequency
        for j=1:size(moviepow,2)
            for k=1:size(moviepow,1)
                minclust = clusti(1) ; maxclust = clusti(end) ; 
                cluststep =  ((moviepow(k,j,maxclust)-moviepow(k,j,minclust)) / length(clusti)) ; 
                if cluststep ~= 0
                    linevals = moviepow(k,j,minclust):cluststep:moviepow(k,j,maxclust) ; 
                    linevals = linevals(1:length(clusti)) ; 
                    newmovie(k,j,clusti) = linevals ; 
                else
                    linevals = ones(1,length(clusti)) * minclust ; 
                    newmovie(k,j,clusti) = linevals ; 
                end            
            end           
        end   
    end
    
    restinds = 100:5:350 ; 
    movieinds = 100:10:635 ; 

    cd(['c:/shared/badger_mri/',subs{sub},'/nii/antsf']) ; 
    
    rest = dir('warped_resting.nii.gz') ; 
    rest = load_untouch_nii(rest(1).name) ; 
    rest = rest.img ; 
   
    movie = load_untouch_nii('../warped_movie.nii.gz') ; 
    movie = movie.img ; 
    
    meanmoviepow = squeeze(mean(newmovie(elecs,:,:))) ; 
    meanrestpow = squeeze(mean(newrest(elecs,:,:))) ; 
    %meanrestpow = eegfiltfft(meanrestpow,1./0.693,.02,1.5) ; 
    for i=1:50 ; disp(i) ;
       voxmoviecorrs(:,:,:,i) = voxcorr(movie(:,:,:,40:end-40),meanmoviepow(i,30:end-50)) ;  
       %voxrestcorrs(:,:,:,i) = voxcorr(rest(:,:,:,40:end-40),meanrestpow(i,30:end-50)) ;  
    end
    
    subvoxmoviecorrs(sub,:,:,:,:) = voxmoviecorrs ; 
    %subvoxrestcorrs(sub,:,:,:,:) = voxrestcorrs ; 

   figure, for i=1:50  ;subplot(5,10,i) ; imagesc(squeeze(mean(voxmoviecorrs(:,:,8:12,i),3)),[-.25,.25]) ; end

end


cd c:/shared/frefs ; sums = load_untouch_nii('sumsumreg.nii.gz') ; figure, 
binanat = sums.img > mean(sums.img(:)) ; 
for i=7:18  ; subplottight(5,6,i) ;
    img = squeeze(mean(mean(subvoxmoviecorrs(:,:,:,i,6:12),5),1)) ; 
    anat = squeeze(sums.img(:,:,i)) ; 
    img = img.*binanat(:,:,i) ; 
    plotoverlayIntensity2D(anat,(abs(img).*15),mat2gray(img),270) ; 
end
