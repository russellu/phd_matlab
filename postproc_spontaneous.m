cd c:/shared/save_paper ; ls 
moviecorrs = load('mean_moviecorrs.mat') ; moviecorrs = moviecorrs.meancorrs ; 
restcorrs = load('mean_restcorrs.mat') ; restcorrs = restcorrs.meancorrs ; 

moviebold = load('movie_sb.mat') ; moviebold = moviebold.allsb ; 
restbold = load('rest_sb.mat') ; restbold = restbold.allsb ; 

moviespecs = load('movie_specs.mat') ; moviespecs = moviespecs.meanspecs ; 
restspecs = load('rest_specs.mat') ; restspecs = restspecs.meanspecs ; 

for i=1:100
    for j=1:40
        [h,p,ci,stats] = ttest(squeeze(restcorrs(:,i,j)).^2,squeeze(moviecorrs(:,i,j)).^2) ; 
        allts(i,j) = stats.tstat ; 
        allps(i,j) = p ; 
    end
end

tmoviecorrs = mean(moviecorrs(:,:,27:34),3) ;
trestcorrs = mean(restcorrs(:,:,27:34),3) ; 

for i=1:size(tmoviecorrs,2)
   [h,p,ci,stats] = ttest(tmoviecorrs(:,i),trestcorrs(:,i)) ;  
   allstats(i) = stats.tstat ; fps(i) = p ; 
end



trtimes = -19:20 ; times = trtimes*0.693 ; 

restcorrs(isnan(restcorrs)) = 0 ; 
subplot(3,2,1) ; imagesc(times,1:100,squeeze(mean(restcorrs,1)),[-.3,.3]) ; axis xy ; vline(0,'k') ; xlabel('time(s)') ; ylabel('frequency(hz)') ; title('REST (eyes closed)') ; 
subplot(3,2,2) ; imagesc(times,1:100,squeeze(mean(moviecorrs,1)),[-.3,.3]) ; axis xy ; vline(0,'k') ;  title('MOVIE (eyes open)') ; 
subplot(3,2,3) ; 
shadedErrorBar([],squeeze(mean(mean(restcorrs(:,8:25,:),1),2)),squeeze(std(mean(restcorrs(:,8:25,:),2),0,1))/sqrt(8),{'r'}) ; hold on ; 
shadedErrorBar([],squeeze(mean(mean(moviecorrs(:,8:25,:),1),2)),squeeze(std(mean(moviecorrs(:,8:25,:),2),0,1))/sqrt(8),{'b'}) ; set(gca,'XTick',1:5:40,'XTickLabel',round(times(1:5:40))) ; 
mrest = squeeze(mean(mean(restcorrs(:,8:25,:),1),2)) ; minrest = find(mrest==min(mrest)) ; vline(minrest,'r') ; text(minrest,0.1,[num2str(times(minrest)),'s']) ; 
mmovie = squeeze(mean(mean(moviecorrs(:,8:25,:),1),2)) ; minmovie = find(mmovie==min(mmovie)) ; vline(minmovie,'b') ; text(minmovie,0.06,[num2str(times(minmovie)),'s']) ; 
ylim([-.3,.2]) ; xlim([1,40]) ; hline(0,'k') ; vline(20,'k') ; xlabel('time(s)') ; ylabel('coupling(r)') ; 
title('alpha/beta (8-25Hz) temporal coupling') ; 

subplot(3,2,4) ; freqz = 40:90 ; 
shadedErrorBar([],squeeze(mean(mean(restcorrs(:,freqz,:),1),2)),squeeze(std(mean(restcorrs(:,freqz,:),2),0,1))/sqrt(8),{'r'}) ; hold on ; 
shadedErrorBar([],squeeze(mean(mean(moviecorrs(:,freqz,:),1),2)),squeeze(std(mean(moviecorrs(:,freqz,:),2),0,1))/sqrt(8),{'b'}) ; set(gca,'XTick',1:5:40,'XTickLabel',round(times(1:5:40))) ; 
%mrest = squeeze(mean(mean(restcorrs(:,freqz,:),1),2)) ; minrest = find(mrest==max(mrest)) ; vline(minrest,'r') ; text(minrest,0.1,[num2str(times(minrest)),'s']) ; 
mmovie = squeeze(mean(mean(moviecorrs(:,freqz,:),1),2)) ; minmovie = find(mmovie==max(mmovie)) ; vline(minmovie,'m') ; text(minmovie,0.09,[num2str(times(minmovie)),'s']) ; 
ylim([-.1,.1]) ; xlim([1,40]) ; hline(0,'k') ; vline(20,'k') ; xlabel('time(s)') ; ylabel('coupling(r)') ; 
title('gamma (40-90Hz) temporal coupling') ; 



subplot(3,2,6) ;
shadedErrorBar([],squeeze(mean(mean(restcorrs(:,:,27:34),1),3)),squeeze(std(mean(restcorrs(:,:,27:34),3),0,1))/sqrt(8),{'r'}) ; hold on ; 
shadedErrorBar([],squeeze(mean(mean(moviecorrs(:,:,27:34),1),3)),squeeze(std(mean(moviecorrs(:,:,27:34),3),0,1))/sqrt(8),{'b'}) ; hline(0,'k')
for i=1:length(fps) ; if fps(i) < 0.01 ; text(i,0.12,'*') ; end ; end ; xlabel('frequency(hz)') ; ylabel('coupling(r)') ; ylim([-.25,.15])
title('frequency coupling at 6-10s time lag, *p<0.01,uncorrected') ; 




imagesc([-.3,.3]) ; colorbar ; 
plot([0,1]) ; hold on ; plot([1,0],'r') ; legend({'movie(EO)','rest(EC)'}) 


smthrest = (imfilter(restspecs,fspecial('gaussian',[1,10],5))) ; smthrest(:,end-10:end) = smthrest(:,end-20:end-10) ; 
smthmovie = (imfilter(moviespecs,fspecial('gaussian',[1,10],5))) ; smthmovie(:,end-10:end) = smthmovie(:,end-20:end-10) ; 


mrest = squeeze(mean(restcorrs(:,:,27:40),3)) ; 
mmovie = squeeze(mean(moviecorrs(:,:,27:40),3)) ;

specdiff = mean(smthrest-smthmovie) ; resdiff = imresize(specdiff(1:425),[1,100]) ; 
resspec = imresize((smthrest + smthmovie)/2,[8,100]) ; resrest = imresize(smthrest,[8,100]) ; resmovie = imresize(smthmovie,[8,100]) ; 
shadedErrorBar([],mean(resrest),std(resrest,0,1)/sqrt(8),{'r'}) ; hold on ; 
shadedErrorBar([],mean(resmovie),std(resmovie,0,1)/sqrt(8),{'b'}) ; hold on ;  ylabel('log power (uv^2)') ; xlabel('frequency (hz)') ; 


smthdiff = resrest-resmovie ; 
coupdiff = mean(mrest)-mean(mmovie) ; subcoupdiff = mrest-mmovie ; 
subplot(2,2,1) ; shadedErrorBar([],mean(subcoupdiff),std(subcoupdiff,0,1)/sqrt(8)) ; hline(0,'k') ; xlabel('frequency(hz)') ; ylabel('coupling difference (rho)') ; 
subplot(2,2,3) ; shadedErrorBar([],mean(smthdiff),std(smthdiff,0,1)/sqrt(8)) ; hline(0,'k') ; xlabel('frequency (hz)') ; ylabel('power difference (log uv^2)') ; ylim([-1,13])
[c,p] = corr(coupdiff(1:100)',resdiff(1:100)') ; 
subplot(1,2,2) ; 
plot(coupdiff(1:100),resdiff(1:100),'o') ; title(['rho=',num2str(c),', rho^2=',num2str(c^2),', p=',num2str(p)]) ; lsline ; 


mboth = (mrest + mmovie) / 2 ; subplot(2,2,1) ; 
shadedErrorBar([],mean(mboth),std(mboth,0,1)/sqrt(8)) ; hline(0,'k') ; ylim([-.2,.1]) ; xlabel('frequency (hz)') ; ylabel('coupling (rho)') ; 
resspec(:,100) = resspec(:,99) ; subplot(2,2,3) ; 
shadedErrorBar([],mean(resspec),std(resspec,0,1)/sqrt(8)) ; ylim([-30,0]) ; ylabel('log power') ; xlabel('frequency (hz)') ; 
subplot(1,2,2) ; plot(mean(resspec),mean(mboth),'o') ; 

mov = load('movie_allmfilts.mat') ; mov = mov.allmfilts ; 
rest = load('rest_allmfilts.mat') ; rest = rest.allmfilts ; 
both = (mov + rest)/2 ; 
subplot(1,2,2) ; plot(mean(resspec),abs(mean(mboth)),'o') ; lsline ; 
[c,p] = corr(mean(resspec)',abs(mean(mboth))') ; 
title(['rho=',num2str(c),', rho^2=',num2str(c^2),', p=',num2str(p)]) ; xlabel('log power') ; ylabel('coupling (abs. rho)') ; 
