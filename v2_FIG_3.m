[mnis,corrs] = getMNI ; 
[meancomps,fmri] = isolateComponent ;times = load('times') ; times = times.times ;
meancomps = meancomps(:,:,1:110,:) ; 
stimtypes = {[1,3,2],[1,5,6]} ;
typenames = {'contrast','randomization'} ;
allnames = {'unperturbed','5%contrast','33%contrast','plaid','10%randomization','60%randomization'} ;
cnames = {{'100%contrast','33%contrast','5%contrast'},{'0%rand','10%rand','60%rand'}} ; cticks = {1:3,1:3} ; 
stimcolors = {[1,0,0],[0,1,0],[0,0,1]} ; 
ctitles = {'contrast%','randomization%'} ;
bandnames = {'alpha(8:15Hz)','beta(15:25Hz','low gamma(40:60Hz)','mid gamma(60:80Hz)','high gamma(80:100Hz)'} ; 
bandnames2 = {'alpha','beta','lowgamma','midgamma','highgamma'} ; 
freqbands = {8:15,15:25,40:60,60:80,80:100} ; 
meanf = squeeze(mean(fmri(:,:,3:8),3))*100 ; 
mmcomps = squeeze(mean(meancomps,1))  ;
stiminds = [2,3,1,5,6] ; 
mtcomps = squeeze(mean(meancomps(:,:,:,times<2 & times >.5),4)) ; 

stimcomps = squeeze(mean(meancomps,1)) ; 
meanf = squeeze(mean(fmri(:,:,4:8),3)) ; 
clear corrmat ; 
for i=1:size(stimcomps,2)
    for j=1:size(stimcomps,3) ;
        e1 = stimcomps(stiminds,i,j) ;
        f1 = mean(meanf(:,stiminds),1) ; 
       [rho,p] = corr([e1,f1']) ; 
       corrmat(i,j) = rho(1,2) ; 
       pmat(i,j) = p(1,2) ; 
    end
end
% get the subject specific scatter plots as well.  yay.
clear subcorrmat ; 
for i=1:size(meancomps,3) ; 
    for j=1:size(meancomps,4) 
        evec = reshape(squeeze(meancomps(:,stiminds,i,j)),[1,size(meancomps,1).*size(stiminds,2)]) ; 
        fvec = reshape(meanf(:,stiminds),[1,size(meancomps,1).*size(stiminds,2)]) ; 
        both = [evec;fvec]' ; 
        [rho,pval] = corr(both) ; 
        allrhos(i,j) = rho(1,2) ; 
        allpvals(i,j) = pval(1,2) ; 
    end
end


stt = find(times>0) ; stt = stt(1) ; endt = find(times>2) ; endt = endt(1) ; mts = min(times):.25:max(times)  ; 
sigspec = flipud((pmat<0.05).*corrmat) ; 
subplot(2,5,1) ; imagesc(sigspec) ;set(gca,'Ytick',1:10:110) ; set(gca,'YTickLabel',110:-10:1) ; vline([stt,endt],'k') ; xlabel('time(s)') ; ylabel('frequency(hz)') ;  set(gca,'XTick',[stt,endt]) ; set(gca,'XTickLabel',[0,2]) ;
%subplot(2,5,2) ; imagesc(sigspec) ;set(gca,'Ytick',1:10:110) ; set(gca,'YTickLabel',110:-10:1) ; vline([stt,endt],'k') ; xlabel('time(s)') ; ylabel('frequency(hz)') ;  set(gca,'XTick',[stt,endt]) ; set(gca,'XTickLabel',[0,2]) ; colorbar
low = 10:30 ; high = 75:110 ; 
lowsubs = squeeze(mean(mean(meancomps(:,:,low,times<2&times>0.5),3),4)) ; highsubs = squeeze(mean(mean(meancomps(:,:,high,times<2&times>0.5),3),4)) ;
mboth = meanf ;
subplot(2,5,3) ; y1 = squeeze(mean(mboth(:,stims),1)) ; x1 = squeeze(mean(lowsubs(:,stims),1)) ; P = polyfit(x1,y1,1); yfit = P(1)*x1+P(2); plot(x1,yfit,'k','LineWidth',2); hold on ; 
errorbarxy(mean(lowsubs(:,stims),1),mean(mboth(:,stims),1),std(lowsubs(:,stims),0,1)./sqrt(22),std(mboth(:,stims),0,1)./sqrt(22),{'o','k','k'}) ; ylabel('BOLD %change in conjunction ROI') ; xlabel('EEG 10-30Hz (db)') ; 
x1 = double(x1) ; y1 = double(y1) ; 
text(x1, y1,labels,'horizontal','left', 'vertical','bottom') ;
[rho,pval] = corr([mean(lowsubs(:,stims),1);mean(mboth(:,stims),1)]') ; 
title(['rho = ',num2str(rho(1,2)),' p = ',num2str(pval(1,2))]) ;
subplot(2,5,2) ; lowsubs = highsubs ;  y1 = squeeze(mean(mboth(:,stims),1)) ; x1 = squeeze(mean(lowsubs(:,stims),1)) ; P = polyfit(x1,y1,1); yfit = P(1)*x1+P(2); plot(x1,yfit,'k','LineWidth',2); hold on ; 
errorbarxy(mean(lowsubs(:,stims),1),mean(mboth(:,stims),1),std(lowsubs(:,stims),0,1)./sqrt(22),std(mboth(:,stims),0,1)./sqrt(22),{'o','k','k'}) ; ylabel('BOLD %change in conjunction ROI') ; xlabel('EEG 75-110Hz (db)') ; 
x1 = double(x1) ; y1 = double(y1) ; 
text(x1, y1,labels,'horizontal','left', 'vertical','bottom') ;
[rho,pval] = corr([mean(lowsubs(:,stims),1);mean(mboth(:,stims),1)]') ; 
title(['rho = ',num2str(rho(1,2)),' p = ',num2str(pval(1,2))]) ;

sigspec = flipud(allrhos.*(allpvals<.05)) ; 
subplot(2,5,6) ; imagesc(sigspec) ;set(gca,'Ytick',1:10:110) ; set(gca,'YTickLabel',110:-10:1) ; vline([stt,endt],'k') ; xlabel('time(s)') ; ylabel('frequency(hz)') ;  set(gca,'XTick',[stt,endt]) ; set(gca,'XTickLabel',[0,2]) ;

evec = reshape(squeeze(mean(mean(meancomps(:,stiminds,high,times<2&times>.5),3),4)),[1,size(meancomps,1).*size(stiminds,2)]) ; 
fvec = reshape(meanf(:,stiminds),[1,size(meancomps,1).*size(stiminds,2)]) ; 
[pi,p,ol] = Shepherd(evec',fvec',10000) ;
subplot(2,5,7) ; x1 = evec(ol==0) ; y1 = fvec(ol==0) ; P = polyfit(x1,y1,1); yfit = P(1)*x1+P(2); plot(x1,yfit,'k','LineWidth',2); hold on ; 
plot(evec(ol==0),fvec(ol==0),'o') ; title(['Shepherds pi = ',num2str(pi),', p = ',num2str(p)]) ;ylabel('BOLD %change') ; xlabel('EEG 75-110Hz (db)') ; 

evec = reshape(squeeze(mean(mean(meancomps(:,stiminds,low,times<2&times>.5),3),4)),[1,size(meancomps,1).*size(stiminds,2)]) ; 
fvec = reshape(meanf(:,stiminds),[1,size(meancomps,1).*size(stiminds,2)]) ; 
[pi,p,ol] = Shepherd(evec',fvec',10000) ;
subplot(2,5,8) ; x1 = evec(ol==0) ; y1 = fvec(ol==0) ; P = polyfit(x1,y1,1); yfit = P(1)*x1+P(2); plot(x1,yfit,'k','LineWidth',2); hold on ; 
plot(evec(ol==0),fvec(ol==0),'o') ; title(['Shepherds pi = ',num2str(pi),', p = ',num2str(p)]) ;ylabel('BOLD %change') ; xlabel('EEG 10-30Hz (db)') ; 

subplot(1,2,1) ; imagesc([-.5,.5]) ; colorbar ; subplot(1,2,2) ; imagesc([-1,1]) ; colorbar ; 

for x=1:110 ; 
    subplot(11,11,x) ; 
for i=1:max(size(stiminds)) ;
   hold on ; 
   plot(mean(meanf(:,stiminds(i)),1),squeeze(mean(mean(mtcomps(:,stiminds(i),x),3),1)),'o','Color',[i/5,0,0],'LineWidth',2) ; 
    
    
    
    
    
end
end












