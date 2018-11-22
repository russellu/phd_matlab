clear all  ; close all
subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie','tegan'} ;
ocomps = {[5,26,32,42],[8,21,23,31],[4,14,15,43],[5,11,23],[5,7,11],[3,11,22],[7,11,15]} ; 
icomps = {[17,21,26,39],[14,19,33],[8,18,24,32],[9,25],[15,17,19],[12,22,28],[7,19]} ; 

for sub=1:length(subs) ; 
    cd(['c:/shared/badger_eeg/',subs{sub},'/outside']) ; ls  ; 
    if sub==1 ; otimes = load('times') ; otimes = otimes.times ; freqs = load('freqs') ; freqs = freqs.freqs ;  end
    bersp = load('bersp.mat') ; bersp = bersp.bersp ; bersp = squeeze(bersp(2,:,1:15,:,:)) ; 
    allbersp(1,sub,:,:,:,:) = bersp ; 
    eego = pop_loadset('eegi.set') ; 
    alleego{sub} = eego ; 
    %figure,
    %for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(bersp(i,:,:,:),2)),[-6,6]) ; axis xy ; title(i) ; end ; suptitle(['subs_',subs{sub},'_outside']) ; 
    subplot(2,7,sub) ; imagesc(squeeze(mean(mean(bersp(ocomps{sub},:,:,:),1),2)),[-6,6]) ; 
    
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    if sub==1 ; itimes = load('times') ; itimes = itimes.times ; end
    bersp = load('bersp.mat') ; bersp = bersp.bersp ; bersp = squeeze(bersp(1,:,1:15,:,:)) ; 
    allbersp(2,sub,:,:,:,:) = bersp ; 
    eegi = pop_loadset('eegi.set') ; 
    alleegi{sub} = eegi ; 
    %figure,
    %for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(bersp(i,:,:,:),2)),[-6,6]) ; axis xy ; title(i) ; end ; suptitle(['subs_',subs{sub},'_inside']) ; 
    subplot(2,7,sub+7) ; imagesc(squeeze(mean(mean(bersp(icomps{sub},1:15,:,:),1),2)),[-6,6]) ; 
end
%{
% process the topoplots
for i=1:length(alleegi)
   figure,
   for j=1:length(icomps{i}) ; 
       subplot(2,5,j) ; topoplot(squeeze(alleegi{i}.icawinv(:,icomps{i}(j))),alleegi{i}.chanlocs) ;       
   end
   for j=1:length(ocomps{i}) ;
       subplot(2,5,j+5) ; topoplot(squeeze(alleego{i}.icawinv(:,ocomps{i}(j))),alleegi{i}.chanlocs) ;       
   end    
end

% find the correlated components
clear corrs ps
for i=1:length(alleegi)
    clear compis compos ; 
   for j=1:length(icomps{i})
      compis(j,:) = alleegi{i}.icawinv(:,icomps{i}(j)) ;  
       
   end
   for j=1:length(ocomps{i})
      compos(j,:) = alleego{i}.icawinv(:,ocomps{i}(j)) ;        
   end
   [corrs{i},ps{i}] = corr(compis',compos') ; 
    
end
% get the indices
figure,
for i=1:length(corrs)
   subplot(2,4,i), imagesc(corrs{i}.^2) ;  colorbar ; 
    
end

ocomps = {[5,26,32,42],[8,21,23,31],[4,14,15,43],[5,11,23],[5,7,11],[3,11,22],[7,11,15]} ; 
icomps = {[17,21,26,39],[14,19,33],[8,18,24,32],[9,25],[15,17,19],[12,22,28],[7,19]} ; 

ininds = {[1,4],[1,2,3],[3,4,2],[1,2],[1,2],[1,2,3],[1,2]} ; 
in2out = {[1,4],[1,4,3],[2,4,1],[1,3],[1,2],[2,3,1],[1,2]} ; 
polarities = {[1,-1],[1,1,1],[-1,1,1],[1,-1],[1,1],[1,1,1],[1,-1]} ; 

for i=1:7
    figure,
    for j=1:length(ininds{i})
       subplot(2,3,j) ; 
       %topoplot(alleegi{i}.icawinv(:,icomps{i}(ininds{i}(j)))*polarities{i}(j),alleegi{i}.chanlocs) ; 
       imagesc(squeeze(mean(allbersp(2,i,icomps{i}(ininds{i}(j)),:,:,:),4))) ; 
       subplot(2,3,j+3) ; 
       imagesc(squeeze(mean(allbersp(1,i,ocomps{i}(in2out{i}(j)),:,:,:),4))) ; 
       %topoplot(alleego{i}.icawinv(:,ocomps{i}(in2out{i}(j))),alleego{i}.chanlocs) ; 
    end
end

%}
% by sorting the top modulation only
% inside
for j=1:7 ; figure,
   in_mij = squeeze(mean(mean(allbersp(2,j,icomps{j},:,:,itimes>0 & itimes<5),4),6)) ; 
   inmod = squeeze(mean(abs(in_mij),2)) ; 
   [~,si] = sort(inmod,'descend') ; 
   for m=1:length(icomps{j})
       subplot(2,5,m) ; 
       imagesc(squeeze(mean(allbersp(2,j,icomps{j}(si(m)),:,:,:),4)),[-5,5]) ; 
       meancomps(1,j,:,:,:) = squeeze(mean(allbersp(2,j,icomps{j}(si(1:2)),:,:,:),3)) ; 
   end
   
   out_mij = squeeze(mean(mean(allbersp(1,j,ocomps{j},:,:,otimes>0 & otimes<3),4),6)) ; 
   outmod = squeeze(mean(abs(out_mij),2)) ; 
   [~,si] = sort(outmod,'descend') ; 
   for m=1:length(ocomps{j})
       subplot(2,5,m+5) ; 
       imagesc(squeeze(mean(allbersp(1,j,ocomps{j}(si(m)),:,:,:),4)),[-5,5]) ; 
       meancomps(2,j,:,:,:) = squeeze(mean(allbersp(1,j,ocomps{j}(si(1:2)),:,:,:),3)) ; 
   end
end
mtcomps(1,:,:,:) = squeeze(mean(meancomps(1,:,:,:,itimes>0 & itimes<3),5)) ; 
mtcomps(2,:,:,:) = squeeze(mean(meancomps(2,:,:,:,otimes>0 & otimes<3),5)) ; 




for i=1:7 ;
   subplot(4,7,i) ; imagesc(itimes,freqs,squeeze(mean(meancomps(1,i,:,:,:),3)),[-5,5]) ; axis xy ; title(['subject ',num2str(i),' inside']) ; vline([0,5],'k') ;
   if i==1 ; xlabel('time(s)') ; ylabel('frequency(hz)') ; end
   subplot(4,7,i+7) ; imagesc(otimes,freqs,squeeze(mean(meancomps(2,i,:,:,:),3)),[-5,5]) ; axis xy ; title(['subject ',num2str(i),' outside']) ; vline([0,3],'k') ; 
   subplot(2,7,i+7) ;  
   errorbar(squeeze(mean(mtcomps(1,i,:,:),3)),squeeze(std(mtcomps(1,i,:,:),0,3))./sqrt(15)) ; hold on ; xlim([0,61]) ; ylim([-7,7]) ; 
   errorbar(squeeze(mean(mtcomps(2,i,:,:),3)),squeeze(std(mtcomps(2,i,:,:),0,3))./sqrt(15),'r') ; xlim([0,61]) ; ylim([-7,7]) ;  hline(0,'k') ; 
   if i==1 ; xlabel('frequency(hz)') ; ylabel('power(db)') ; end
   set(gca,'XTick',1:10:60,'XTickLabel',round(freqs(1:10:60))) ; title(['subject ',num2str(i),' out(red) vs in(blue)']) ;
end

for i=1:60  
[h,p,ci,stats] = ttest(squeeze(mean(mtcomps(1,:,:,i),3)),squeeze(mean(mtcomps(2,:,:,i),3))) ;
allstats(i) = stats.tstat ; 
allps(i) = p ; 
end
% grand average
errorbar(squeeze(mean(mean(mtcomps(1,:,:,:),2),3)),squeeze(std(mean(mtcomps(1,:,:,:),3),0,2))./sqrt(7)) ; hold on ; 
errorbar(squeeze(mean(mean(mtcomps(2,:,:,:),2),3)),squeeze(std(mean(mtcomps(2,:,:,:),3),0,2))./sqrt(7),'r') ; 
for i=1:length(allps) ; if allps(i) < 0.05 ; text(i,3,'*') ; end ; end
xlim([0,61]) ; ylim([-5.5,3.5]) ;  hline(0,'k') ; 
xlabel('frequency(hz)') ; ylabel('power(db)') ; 
set(gca,'XTick',1:10:60,'XTickLabel',round(freqs(1:10:60))) ; title('grand average inside vs outside, *p<0.05') ; 
legend({'inside','outside'}) ; 

% correlations : alpha/beta, and gamma. show scatter plots
figure,
freqintrs = freqs>35 & freqs<80; 
errorbarxy(squeeze(mean(mean(mtcomps(1,:,:,freqintrs),3),4)),squeeze(mean(mean(mtcomps(2,:,:,freqintrs),3),4)),...
squeeze(std(mean(mtcomps(1,:,:,freqintrs),4),0,3))./sqrt(7),squeeze(std(mean(mtcomps(2,:,:,freqintrs),4),0,3))./sqrt(7),{'o','k','k'}) ; hold on ; 
[c,p] = corr([squeeze(mean(mean(mtcomps(1,:,:,freqintrs),3),4));squeeze(mean(mean(mtcomps(2,:,:,freqintrs),3),4))]') ;
title(['r=',num2str(c(1,2)),' p=',num2str(p(1,2))]) ; xlim([0,3]) ; ylim([0,3]) ; suptitle('gamma (35-80Hz) power inside vs outside') ; xlabel('35-80Hz power (inside') ; ylabel('35-80Hz power (outside)') ; 
figure,
freqintrs = freqs>10 & freqs<25; 
errorbarxy(squeeze(mean(mean(mtcomps(1,:,:,freqintrs),3),4)),squeeze(mean(mean(mtcomps(2,:,:,freqintrs),3),4)),...
squeeze(std(mean(mtcomps(1,:,:,freqintrs),4),0,3))./sqrt(7),squeeze(std(mean(mtcomps(2,:,:,freqintrs),4),0,3))./sqrt(7),{'o','k','k'}) ; hold on ; 
[c,p] = corr([squeeze(mean(mean(mtcomps(1,:,:,freqintrs),3),4));squeeze(mean(mean(mtcomps(2,:,:,freqintrs),3),4))]') ;
title(['r=',num2str(c(1,2)),' p=',num2str(p(1,2))]) ; xlim([-5,1]) ; ylim([-5,1]) ; suptitle('alpha/beta (10-25Hz) power inside vs outside') ; xlabel('10-25Hz power (inside') ; ylabel('10-25Hz power (outside)') ; 

% single trials: two separate subjects, and all subjects t-tests.
sub=2 ; 
figure,
for i=1:15 ; 
    subplot(2,15,i) ; imagesc(itimes,freqs,imfilter(squeeze(meancomps(1,sub,i,:,:)),fspecial('gaussian',[3,7],3)),[-10,10]) ; axis xy ; if i==1 ; xlabel('time(s)') ; ylabel('frequency(hz)') ; end ; vline([0,5],'k') ; 
    subplot(2,15,i+15) ; imagesc(otimes,freqs,imfilter(squeeze(meancomps(2,sub,i,:,:)),fspecial('gaussian',[3,7],3)),[-10,10]) ; axis xy ; vline([0,3],'k')  ;
end

for i=1:7 ; 
   subplot(2,7,i) ; plot(squeeze(mtcomps(1,i,:,:))','b') ; ylim([-15,10]) ; hline(0,'k') ; xlim([0,61]) ; title(['subject ',num2str(i),' out(red) vs in(blue)']) ;
   if i==1 ; xlabel('frequency(hz)') ; ylabel('power(db)') ; end
   set(gca,'XTick',1:10:60,'XTickLabel',round(freqs(1:10:60))) 
   subplot(2,7,i+7) ; plot(squeeze(mtcomps(2,i,:,:))','r') ;  ylim([-15,10]) ; hline(0,'k') ; xlim([0,61]) ; 
   set(gca,'XTick',1:10:60,'XTickLabel',round(freqs(1:10:60))) 
end

clear ps tstats 
for i=1:7
    for j=1:60
        [h,p,ci,stat] = ttest2(squeeze(mtcomps(1,i,:,j)),squeeze(mtcomps(2,i,:,j))) ; 
        ps(i,j) = p ; tstats(i,j) = stat.tstat ; 
    end
end


errorbar(mean(tstats,1),std(tstats,0,1)./sqrt(7),'k') ; hline(0,'r') ; 
set(gca,'XTick',1:10:60,'XTickLabel',round(freqs(1:10:60))) ; ylabel('T value (inside-outside)') ; xlabel('frequency(hz)') ; xlim([0,61]) ; 
suptitle('frequency specific T values for inside-outside') ; 

%stdcomps = squeeze(std(mtcomps,0,3)) ; 
%subplot(1,2,1) ; imagesc(squeeze(stdcomps(1,:,:)),[0,5]) ; 
%subplot(1,2,2) ; imagesc(squeeze(stdcomps(2,:,:)),[0,5]) ; 
errorbar(squeeze(mean(stdcomps(1,:,:),2)),squeeze(std(stdcomps(1,:,:),0,2))./sqrt(7)) ; hold on  ;
errorbar(squeeze(mean(stdcomps(2,:,:),2)),squeeze(std(stdcomps(2,:,:),0,2))./sqrt(7),'r') ;
set(gca,'XTick',1:10:60,'XTickLabel',round(freqs(1:10:60))) ; ylabel('inter-trial power std (db)') ; xlabel('frequency(hz)') ; xlim([0,61]) ; 
legend({'inside','outside'}) ; suptitle('') ; 
suptitle('mean inter-trial standard deviation across subjects') ; 


