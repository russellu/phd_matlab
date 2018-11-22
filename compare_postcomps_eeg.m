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
    
    for i=1:size(restpow,1)
        for j=1:size(restpow,2)
             for k=1:4
                   shiftcount = 1 ;         
                for shift = -20:20
                    
                    if k==4
                        [c,p] = corr([circshift(squeeze(restpow(i,j,40:end-40)),shift),mean(segmix{2}(40:end-40,compinds{sub}(k:end)),2)]) ;
                        restcorrs(sub,i,j,shiftcount,k) = c(1,2) ; restps(sub,i,j,k) = p(1,2) ; 
                        moviecorrs(sub,i,j,shiftcount,k) = corr2(circshift(squeeze(moviepow(i,j,40:end-40)),shift),mean(segmix{1}(40:end-40,compinds{sub}(k:end)),2)) ; 
                        shiftcount = shiftcount + 1 ;          
                    else             
                        [c,p] = corr([circshift(squeeze(restpow(i,j,40:end-40)),shift),segmix{2}(40:end-40,compinds{sub}(k))]) ;
                        restcorrs(sub,i,j,shiftcount,k) = c(1,2) ; restps(sub,i,j,k) = p(1,2) ; 
                        moviecorrs(sub,i,j,shiftcount,k) = corr2(circshift(squeeze(moviepow(i,j,40:end-40)),shift),segmix{1}(40:end-40,compinds{sub}(k))) ; 
                        shiftcount = shiftcount + 1 ;     
                    end
                end
            end
        end
    end
 
   
    rests=dir('*rest*Pulse*vhdr') ; 
 
    merged = pop_loadset('allstim_broad_merged.set') ;
    lats = {merged.urevent.latency} ; 
    types = {merged.urevent.type} ; 
    r128s = find(strcmp('R128',types)) ; 
    rlats = cell2mat(lats(r128s)) ; 
    start = rlats(1)+merged.srate ; stop = rlats(length(rlats))-merged.srate ; 
    merged = pop_select(merged,'point',[start,stop]) ; 
    tp(merged) ; suptitle(subs{sub}) ; 
    

    

    
end
freqs = 1:2:100 ; times = (-20:1:20) *0.693 ; 

subplot(1,2,1) ; 
imagesc(-20:20,1:2:100,squeeze(mean(mean(mean(restcorrs(:,elecs,:,:,:),1),2),5)),[-.2,.2]) ; vline(0,'k') ; axis xy ;
subplot(1,2,2) ;
imagesc(-20:20,1:2:100,squeeze(mean(mean(mean(moviecorrs(:,elecs,:,:,:),1),2),5)),[-.2,.2]) ; vline(0,'k') ; axis xy ; 

for i=1:7 ; subplot(2,7,i) ; imagesc(times,freqs,squeeze(mean(mean(restcorrs(i,elecs,:,:,[1,2,4]),2),5)),[-.6,.6]) ; axis xy ; vline(0,'k') ; if i==1 ;  xlabel('time(s)')  ;ylabel('freq(hz)') ; end ; end
for i=1:7 ; subplot(2,7,i+7) ; imagesc(times,freqs,squeeze(mean(mean(restcorrs(i,elecs,:,:,[3]),2),5)),[-.6,.6]) ; axis xy ; vline(0,'k') ; if i==1 ;  xlabel('time(s)')  ;ylabel('freq(hz)') ; end ; end

subplot(1,2,1) ; imagesc(times,freqs,squeeze(mean(mean(mean(restcorrs(:,elecs,:,:,[1,2,4]),1),2),5)),[-.25,.15]) ; axis xy ;xlabel('time(s)')  ;ylabel('freq(hz)') ; vline(0,'k') ; title('rest'); 
subplot(1,2,2) ; imagesc(times,freqs,squeeze(mean(mean(mean(moviecorrs(:,elecs,:,:,[1,2,4]),1),2),5)),[-.25,.15]) ; axis xy ; vline(0,'k') ; title('movie') ; 
mrest = squeeze(mean(mean(mean(restcorrs(:,elecs,:,:,[1,2,4]),1),2),5)) ; 
mmovie = squeeze(mean(mean(mean(moviecorrs(:,elecs,:,:,[1,2,4]),1),2),5)) ; 

figure,subplot(1,2,1) ; imagesc([-.25,.15]) ;colorbar ; subplot(1,2,2) ; imagesc([-.6,.6]) ; colorbar

elecmovies = squeeze(mean(mean(moviecorrs(:,elecs,:,:,[1,2,4]),2),5)) ; 
elecrests = squeeze(mean(mean(restcorrs(:,elecs,:,:,[1,2,4]),2),5)) ; 
colors = {'r','k','b','g','c','m','y'} ; 
cfreqs= 1:20; ctimes = 20:41 ; 
for i=1:7
     subplot(1,7,i) ;hold on ; 
    plot(reshape(squeeze(elecrests(i,cfreqs,ctimes)),[1,numel(elecrests(i,cfreqs,ctimes))]),reshape(squeeze(elecmovies(i,cfreqs,ctimes)),[1,numel(elecmovies(i,cfreqs,ctimes))]),['.']) ; 
    title(corr2(reshape(squeeze(elecrests(i,cfreqs,ctimes)),[1,numel(elecrests(i,cfreqs,ctimes))]),reshape(squeeze(elecmovies(i,cfreqs,ctimes)),[1,numel(elecmovies(i,cfreqs,ctimes))]))) ; 
    lsline ;if i==1 ; xlabel('rest') ; ylabel('movie') ; end
end

plot(reshape(mrest(cfreqs,ctimes),[1,numel(mrest(cfreqs,ctimes))]),reshape(mmovie(cfreqs,ctimes),[1,numel(mmovie(cfreqs,ctimes))]),'.') ; lsline 
title(corr2(reshape(mrest(cfreqs,ctimes),[1,numel(mrest(cfreqs,ctimes))]),reshape(mmovie(cfreqs,ctimes),[1,numel(mmovie(cfreqs,ctimes))]))) ; xlabel('rest') ; ylabel('movie') ; 

subplot(1,2,1) ; 
bar([mean(mean(elecrests(:,cfreqs,ctimes).^2,2),3),mean(mean(elecmovies(:,cfreqs,ctimes).^2,2),3)]) ; xlabel('subject') ; ylabel('coupling (mean r^2)') ; legend({'rest','movie'}) ; 
bar(mean(mean(meanfs(:,:,elecs,f>8 & f<12),3),4)) ;
meanpows = mean(mean(meanfs(:,:,elecs,f>8 & f<12),3),4) ;
meancoupling = [mean(mean(elecrests(:,cfreqs,ctimes).^2,2),3),mean(mean(elecmovies(:,cfreqs,ctimes).^2,2),3)] ; 
plot(mean(meancoupling,2),mean(meanpows,2),'o') ; lsline ; [c,p] = corr(mean(meancoupling,2),mean(meanpows,2)) ; title(['corr=',num2str(c),' p=',num2str(p)]) ; xlabel('mean coupling') ; ylabel('mean alpha power') ; 
cd('c:\shared\badger_eeg\alex') ; alldiffs=  load('alldiffs') ; alldiffs = alldiffs.alldiffs ; 
% correlate power with distance and topoplots
for i=1:65 ; corrs(i) = corr2(squeeze(mean(alldiffs(:,i,:),3)),mean(meancoupling,2)) ; end



% get the BOLD power spectrum:
for i=1:length(allsegmix)
    segmovie = allsegmix{i}{1} ;        
    moviecomp = squeeze(mean(segmovie(40:end-40,compinds{i}([1:2,4:end])),2)) ; 
    segrest = allsegmix{i}{2} ;        
    restcomp = squeeze(mean(segrest(40:end-40,compinds{i}([1:2,4:end])),2)) ;    
    
    [restspec(i,:),f] = spectopo(restcomp',0,1./0.693,'plot','off') ; 
    [moviespec(i,:),~] = spectopo(moviecomp',0,1./0.693,'plot','off') ; 

    
end

finds = 5:100 ; 
errorbar(mean(restspec(:,finds),1),std(restspec(:,finds),0,1)./sqrt(7)) ; hold on ; errorbar(mean(moviespec(:,finds),1),std(moviespec(:,finds),0,1)./sqrt(7),'r') ;
set(gca,'XTick',1:10:length(finds),'XTickLabel',round(f(finds(1:10:end))*100)./100) ; xlabel('frequency(hz)') ; ylabel('log power') ; 
suptitle('ALFF range power for movie and rest, *p<0.1') ; legend({'rest','movie'}) ; 
for i=1:257
   [h,p,ci,stats] = ttest(restspec(:,i),moviespec(:,i)) ;  
   ts(i) = stats.tstat ; 
   ps(i) = p ; 
end
pvals = ps(finds) ; text(find(pvals<.1),ones(1,length(find(pvals<.1)))*21,'*') ; 



    %{
freqs = 1:2:100 ; times = (-20:1:20) *0.693 ; 
for s=1:7
    figure,
for i=1:4 ; subplot(2,4,i) ; imagesc(times,freqs,squeeze(mean(mean(restcorrs(s,elecs,:,:,i),1),2)),[-.7,.7]) ; vline(0,'k') ; axis xy ; end ; 
for i=1:4 ; subplot(2,4,i+4) ; imagesc(times,freqs,squeeze(mean(mean(moviecorrs(s,elecs,:,:,i),1),2)),[-.7,.7]) ; vline(0,'k') ; axis xy ; end ; 
suptitle(subs{s}) ; 
end

figure,
for i=1:4 ; subplot(2,4,i) ; imagesc(times,freqs,squeeze(mean(mean(restcorrs(:,elecs,:,:,i),1),2)),[-.35,.35]) ; vline(0,'k') ; axis xy ; end ; 
for i=1:4 ; subplot(2,4,i+4) ; imagesc(times,freqs,squeeze(mean(mean(moviecorrs(:,elecs,:,:,i),1),2)),[-.35,.35]) ; vline(0,'k') ; axis xy ; end ; 
suptitle('grand avg') ; 

subplot(2,2,1) ; 
imagesc(times,freqs,squeeze(mean(mean(mean(restcorrs(:,elecs,:,:,[1,3,4]),1),2),5)),[-.2,.2]) ; axis xy ; vline(0,'k') ; colorbar ; xlabel('time lag(s)') ; ylabel('frequency(hz)') ; 
subplot(2,2,2) ; 
imagesc(times,freqs,squeeze(mean(mean(mean(moviecorrs(:,elecs,:,:,[1,3,4]),1),2),5)),[-.2,.2]) ; axis xy ; vline(0,'k') ; colorbar ; 
mrest = squeeze(mean(mean(mean(restcorrs(:,elecs,:,:,[1,3,4]),1),2),5)) ; 
mmovie = squeeze(mean(mean(mean(moviecorrs(:,elecs,:,:,[1,3,4]),1),2),5)) ;
subplot(2,1,2) ; plot(reshape(mmovie,[1,numel(mmovie)]),reshape(mrest,[1,numel(mrest)]),'.') ; xlabel('movie kernel') ; ylabel('rest kernel') ; 
[c,p] = corr([reshape(mmovie,[1,numel(mmovie)]);reshape(mrest,[1,numel(mrest)])]') ; title(['rho=',num2str(c(1,2)),' p=',num2str(p(1,2))]) ; lsline


melecrest = squeeze(mean(mean(restcorrs(:,:,:,:,[1,3,4]),1),5)) ; 

subplot(1,2,1) ; 
topoplot(squeeze(mean(mean(melecrest(:,7:15,20:25),2),3)),merged.chanlocs,'maplimits',[-.2,.2]) ;
subplot(1,2,2) ; 
topoplot(squeeze(mean(mean(melecrest(:,5,32),2),3)),merged.chanlocs,'maplimits',[-.35,.35]) ;
%}





