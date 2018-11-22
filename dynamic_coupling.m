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
      %{  
    % for rest
    resttrigs = find(strcmp('R128',{restmotion.urevent.type})) ; 
    restlats = cell2mat({restmotion.urevent.latency}) ; 
    restpnts =  restlats(resttrigs) ; 
    restremoveTR = zeros(1,length(restpnts)) ; 
    for i=1:length(restpnts)
        if intersect(restpnts(i),restmotion.badts) 
            restremoveTR(i) = 1 ; 
        end
    end
    
    % for movie
    movietrigs = find(strcmp('R128',{moviemotion.urevent.type})) ; 
    movielats = cell2mat({moviemotion.urevent.latency}) ; 
    moviepnts =  movielats(movietrigs) ; 
    movieremoveTR = zeros(1,length(moviepnts)) ; 
    for i=1:length(moviepnts)
        if intersect(moviepnts(i),moviemotion.badts) 
            movieremoveTR(i) = 1 ; 
        end
    end
    figure,
    subplot(1,3,1) ; imagesc(squeeze(mean(restpow(elecs,20:end,:)))) ; vline(find(restremoveTR==1)) ;
    
    restremoveTR = imdilate(restremoveTR,strel(ones(1,3))) ;    
    movieremoveTR = imdilate(movieremoveTR,strel(ones(1,3))) ; 
    restpow(:,:,restremoveTR==1) = repmat(mean(restpow(:,:,40:end-40),3),[1,1,length(find(restremoveTR==1))]) ; 
    moviepow(:,:,movieremoveTR==1) = repmat(mean(moviepow(:,:,40:end-40),3),[1,1,length(find(movieremoveTR==1))]) ;
 
    allsegmix{sub}{1}(movieremoveTR==1,:) =  repmat(mean(allsegmix{sub}{1}(40:end-40,:),1),[length(find(movieremoveTR==1)),1]) ; 
    allsegmix{sub}{2}(restremoveTR==1,:) = repmat(mean(allsegmix{sub}{1}(40:end-40,:),1),[length(find(restremoveTR==1)),1]) ; 
    %}
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

    for i=1:size(restpow,1)
        for j=1:size(restpow,2)
            for ind=1:length(restinds)
                for k=1:4           
                    if k==4
                        restcorrs(sub,ind,i,j,k,:) = flipud(xcorr(squeeze(newrest(i,j,restinds(ind)-50:restinds(ind)+50)),mean(allsegmix{sub}{2}(restinds(ind)-50:restinds(ind)+50,compinds{sub}(4:end)),2),20,'coeff')) ; 
                    else
                        restcorrs(sub,ind,i,j,k,:) = flipud(xcorr(squeeze(newrest(i,j,restinds(ind)-50:restinds(ind)+50)),allsegmix{sub}{2}(restinds(ind)-50:restinds(ind)+50,compinds{sub}(k)),20,'coeff')) ; 
                    end
                end
            end
        end
    end
    
    for i=1:size(moviepow,1)
        for j=1:size(moviepow,2)
            for ind=1:length(movieinds)
                for k=1:4           
                    if k==4
                        moviecorrs(sub,ind,i,j,k,:) = flipud(xcorr(squeeze(newmovie(i,j,movieinds(ind)-50:movieinds(ind)+50)),mean(allsegmix{sub}{1}(movieinds(ind)-50:movieinds(ind)+50,compinds{sub}(4:end)),2),20,'coeff')) ; 
                    else
                        moviecorrs(sub,ind,i,j,k,:) = flipud(xcorr(squeeze(newmovie(i,j,movieinds(ind)-50:movieinds(ind)+50)),allsegmix{sub}{1}(movieinds(ind)-50:movieinds(ind)+50,compinds{sub}(k)),20,'coeff')) ; 
                    end
                end
            end
        end
    end

end


for s=1:7 ; figure,for i=1:51 ; subplot(6,10,i) ; imagesc(squeeze(mean(restcorrs(s,i,elecs,:,4,:),3))) ; title(i) ; end ; end

freqs = 1:2:100 ; times = (-20:1:20) *0.693 ; 
badrestkerns = {[4:8],[4:15,26:28],[1,3,17,18],[2:6,18,19],[45:50,24:29],[29:33,49:51,7,8,13,14],[11:16,20:24]} ;
badmoviekerns = {[16,21,24,26,32,49,50,51],[1:3,8:16,26,33,35],[2:4,17:18,27:3142:47],[1:5,19:24,42:44,47:49],[1:3,17,21:24,33:36,46],[1:3,12:13,25:34],[27:32,49:51]} ; 
for i=1:size(restcorrs,1) ; subplot(2,9,i) ; 
    goodkerns = zeros(1,size(restcorrs,2)) ; 
    goodkerns(badrestkerns{i}) = 1 ; goodkerns = find(goodkerns==0) ; 
    imagesc(times,freqs,squeeze(mean(mean(restcorrs(i,goodkerns,elecs,:,4,:),2),3)),[-.25,.25]) ; axis xy ; 
    meanrests(i,:,:) = squeeze(mean(mean(restcorrs(i,goodkerns,elecs,:,4,:),2),3)) ; if i==1 ; xlabel('time lag(s)') ; ylabel('frequency(hz)') ; end ; title(['subject ',num2str(i)]) ; vline(0,'k') ; 
end
subplot(2,9,9) ; imagesc(times,freqs,squeeze(mean(meanrests,1)),[-.1,.1]) ; axis xy ; title('grand average (rest EC)') ; vline(0,'k') ; 

for i=1:size(moviecorrs,1) ; subplot(2,9,i+9) ; 
    goodkerns = zeros(1,size(moviecorrs,2)) ; 
    goodkerns(badmoviekerns{i}) = 1 ; goodkerns = find(goodkerns==0) ; 
    imagesc(times,freqs,squeeze(mean(mean(moviecorrs(i,goodkerns,elecs,:,4,:),2),3)),[-.25,.25]) ; axis xy ; 
    meanmovies(i,:,:) = squeeze(mean(mean(moviecorrs(i,goodkerns,elecs,:,4,:),2),3)) ; if i==1 ; xlabel('time lag(s)') ; ylabel('frequency(hz)') ; end ; title(['subject ',num2str(i)]) ; vline(0,'k') ; 
end
subplot(2,9,18) ; imagesc(times,freqs,squeeze(mean(meanmovies,1)),[-.1,.1]) ; axis xy ; title('grand average (movie EO)') ; vline(0,'k') ; 



hz = 4:15 ; 
subplot(1,2,1) ; 
shadedErrorBar([],squeeze(mean(mean(meanrests(:,hz,:),1),2)),squeeze(std(mean(meanrests(:,hz,:),2),0,1))./sqrt(7),'k') ; ylim([-.25,.15]) ; xlim([0,42]) ; hold on ; 
%hz = 20:50 ; shadedErrorBar([],squeeze(mean(mean(meanrests(:,hz,:),1),2)),squeeze(std(mean(meanrests(:,hz,:),2),0,1))./sqrt(7),'k') ; ylim([-.1,.1]) ; xlim([0,42]) ; hline(0,'k') ; 
title('rest coupling') ; ylabel('correlation (r)') ; set(gca,'XTick',1:3:length(times),'XTickLabel',round(times(1:3:end))) ;xlabel('time lag(s)') ; vline(20.5,'r') ; 

%subplot(1,2,2) ;
hz = 4:15 ; 
shadedErrorBar([],squeeze(mean(mean(meanmovies(:,hz,:),1),2)),squeeze(std(mean(meanmovies(:,hz,:),2),0,1))./sqrt(7),'m') ; ylim([-.3,.15]) ; xlim([0,42]) ; hold on ; 
%hz = 20:50 ; shadedErrorBar([],squeeze(mean(mean(meanmovies(:,hz,:),1),2)),squeeze(std(mean(meanmovies(:,hz,:),2),0,1))./sqrt(7),'m') ; 
ylim([-.15,.1]) ; hline(0,'k') ;  
%title('movie coupling') ;
ylabel('correlation (r)') ; set(gca,'XTick',1:3:length(times),'XTickLabel',round(times(1:3:end))) ;xlabel('time lag(s)') ; vline(20.5,'r') ; 

subplot(1,2,2) ; 
shadedErrorBar([],squeeze(mean(mean(meanrests(:,:,31:35),1),3)),squeeze(std(mean(meanrests(:,:,31:35),3),0,1))./sqrt(7),'k') ; hold on ; 
shadedErrorBar([],squeeze(mean(mean(meanmovies(:,:,31:35),1),3)),squeeze(std(mean(meanmovies(:,:,31:35),3),0,1))./sqrt(7),'m') ; hline(0,'k') ; xlim([1,50]) ; 
set(gca,'XTick',1:5:length(freqs),'XTickLabel',round(freqs(1:5:end))) ; xlabel('frequency(hz)') ; ylabel('correlation(r)') ; %legend({'rest','movie'}) ; 
clear h p ci stats
for i=1:50
    [h,p,ci,stats] = ttest(squeeze(mean(meanrests(:,i,31:35),3))',squeeze(mean(meanmovies(:,i,31:35),3))') ;
    pvals(i) = p ; tvals(i) = stats.tstat ; 
end
for i=1:length(pvals) ; if pvals(i)<.05 ; text(i,.05,'*') ; end ; end ; title('* p<0.05, uncorrected') ; 

% get the peaks (rest)
hz = 4:10 ; 
peakinds = find(times>=2 & times<=10) ; 
restcurves = (squeeze(mean(mean(mean(restcorrs(:,:,elecs,hz,4,peakinds),2),3),4))) ;
drest = diff(restcurves,1,2) ; 
clear restpeaks 
for i=1:7
    for j=1:size(drest,2)-1
        if drest(i,j) < 0 && drest(i,j+1) >= 0
            restpeaks(i) = j ;  
        end
    end
end

% get the peaks (movie)
peakinds = find(times>=2 & times<=10) ; peaktimes = times(peakinds) ; 
moviecurves = (squeeze(mean(mean(mean(moviecorrs(:,:,elecs,hz,4,peakinds),2),3),4))) ;
dmovie = diff(moviecurves,1,2) ; 
clear moviepeaks 
for i=1:7
    for j=1:size(dmovie,2)-1
        if dmovie(i,j) < 0 && dmovie(i,j+1) >= 0
            moviepeaks(i) = j ;  
        end
    end
end

bothpeaks = ([restpeaks;moviepeaks] + 3) .*0.693 ; 
[h,p,ci,stats] = ttest(bothpeaks(1,:),bothpeaks(2,:)) ; 
barwitherr(std(bothpeaks,0,2)./sqrt(7),mean(bothpeaks,2)) ; title(['p=',num2str(p),' peak(rest) > peak(movie)']) ; set(gca,'XTickLabel',{'rest','movie'}) ; 
ylabel('lag peak (s)') ; 





subplot(1,3,1) ; imagesc(rands(1:200,1:200,1)) ; colormap gray ; set(gca,'XTickLabel',[],'YTickLabel',[])  ; title('0%SR') ; 
subplot(1,3,2) ; imagesc(rands(1:200,1:200,3)) ; colormap gray ; set(gca,'XTickLabel',[],'YTickLabel',[])  ;  title('10%SR') ; 
subplot(1,3,3) ; imagesc(rand(200,200)) ; colormap gray ; set(gca,'XTickLabel',[],'YTickLabel',[])  ;  title('100%SR') ; 

subplot(1,2,1) ; plot(1,'k') ; hold on ; plot(1,'m') ; legend({'rest (EC)','movie (EO)'}) ; 
subplot(1,2,2) ; plot(1,'r') ; hold on ; plot(1,'b') ; legend({'0%SR and 10%SR','100%SR'}) ; 





