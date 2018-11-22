clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie','tegan'} ; 
freqranges = {[48,60],[45,60],[60,80],[50,70],[50,70],[50,70],[56,70],[55,70]} ; 
lowcomps = {[15,25],[35,30,18],[25,17],[35,20,7],[16,15],[14,33],[25,13,35],[37,14,33]} ;
highcomps = {[26,25,21],[44,37],[39,36],[19,38,50],[20,37],[34,19],[54,24,50],[43,49,46]} ; 

lowcompso = {[25,6,48],[26,8,18],[15,16,5],[12,5],[9,7],[25,10,3],[13,7],[14,12,11]} ;
highcompso = {[20,17,8],[43,25,50],[48,51],[29,14,24],[63,27],[39,30],[32,25],[51,50,49]} ; 

for s=1:length(subs)  ;
    cd(['c:/shared/badger_eeg/',subs{s}]) ;  
    outfiles = dir('*gamma*Pulse*set') ; 
    %disp(outfiles.name) ; 
    EEG = pop_loadset(outfiles(1).name) ; EEG2 = pop_loadset(outfiles(1).name) ; EEG = pop_mergeset(EEG,EEG2) ; 
    EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
    EEG = pop_resample(EEG,250) ; 
    filtlow = EEG ; filtlow.data = eegfiltfft(EEG.data,EEG.srate,1,128) ; 
    filthigh = EEG ; filthigh.data = eegfiltfft(EEG.data,EEG.srate,freqranges{s}(1),freqranges{s}(2)) ; 

    low = pop_loadset('epochlow.set') ; 
    high = pop_loadset('epochhigh.set') ; 
    
    applow = ica_applyweights(EEG,low) ; 
    epochlow1 = pop_epoch(applow,{'S  1'},[-2,7]) ; 
    epochlow2 = pop_epoch(applow,{'S  3'},[-2,7]) ; 
    
    apphigh = ica_applyweights(EEG,high) ; 
    epochhigh1 = pop_epoch(apphigh,{'S  1'},[-2,7]) ; 
    epochhigh2 = pop_epoch(apphigh,{'S  3'},[-2,7]) ; 

    pop_saveset(epochlow1,'epochlow1') ; pop_saveset(epochlow2,'epochlow2') ; 
    pop_saveset(epochhigh1,'epochhigh1') ; pop_saveset(epochhigh2,'epochhigh2') ; 

    
    for i=1:64 ; 
       [ersplow1(i,:,:),~,~,~,~,~,~] = newtimef(squeeze(epochlow1.icaact(i,:,:)),epochlow1.pnts,[epochlow1.xmin,epochlow1.xmax],epochlow1.srate,0,...
           'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',0,'timesout',100,'winsize',64) ;  
       [ersplow2(i,:,:),~,~,~,~,~,~] = newtimef(squeeze(epochlow2.icaact(i,:,:)),epochlow1.pnts,[epochlow1.xmin,epochlow1.xmax],epochlow1.srate,0,...
           'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',0,'timesout',100,'winsize',64) ;  
       [ersphigh1(i,:,:),~,~,~,~,~,~] = newtimef(squeeze(epochhigh1.icaact(i,:,:)),epochhigh1.pnts,[epochhigh1.xmin,epochhigh1.xmax],epochhigh1.srate,0,...
           'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',0,'timesout',100,'winsize',64) ;  
       [ersphigh2(i,:,:),~,~,~,~,~,~] = newtimef(squeeze(epochhigh2.icaact(i,:,:)),epochhigh1.pnts,[epochhigh1.xmin,epochhigh1.xmax],epochhigh1.srate,0,...
           'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',0,'timesout',100,'winsize',64) ;  
    end
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersplow2(i,:,:)),[-10,10]) ; title(i) ; end ; suptitle([subs{s},', low']) ; 
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersphigh2(i,:,:)),[-10,10]) ; title(i) ; end ; suptitle([subs{s},', high']) ; 

end

clear ersplow1 ersplow2 ersphigh1 ersphigh2
for s=1:length(subs)  ;
    cd(['c:/shared/badger_eeg/',subs{s}]) ;  
    low1 = pop_loadset('epochlow1.set') ;     low2 = pop_loadset('epochlow2.set') ; 

    high1 = pop_loadset('epochhigh1.set') ;     high2 = pop_loadset('epochhigh2.set') ; 

    clear ersplow ersphigh
    % for the low frequency filtered data
    for i=1:length(lowcomps{s}) 
        for j=1:15
           [ersplow1(i,j,:,:)] = newtimef(squeeze(low1.icaact(lowcomps{s}(i),:,j)),low.pnts,[low.xmin,low.xmax],low.srate,0,...
               'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',NaN,'timesout',100,'winsize',64) ;   
           [ersplow2(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(low2.icaact(lowcomps{s}(i),:,j)),low.pnts,[low.xmin,low.xmax],low.srate,0,...
               'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',NaN,'timesout',100,'winsize',64) ;   
        end
    end
    % for the high frequency filtered data
    for i=1:length(highcomps{s}) ;
        for j=1:15
           [ersphigh1(i,j,:,:)] = newtimef(squeeze(high1.icaact(highcomps{s}(i),:,j)),high.pnts,[high.xmin,high.xmax],high.srate,0,...
               'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',NaN,'timesout',100,'winsize',64) ;  
           [ersphigh2(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(high2.icaact(highcomps{s}(i),:,j)),high.pnts,[high.xmin,high.xmax],high.srate,0,...
               'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',NaN,'timesout',100,'winsize',64) ;  
        end
    end
    itimes = times ; 
    bersplow1 = ersplow1 - repmat(mean(ersplow1(:,:,:,times<0),4),[1,1,1,100]) ; 
    bersphigh1 = ersphigh1 - repmat(mean(ersphigh1(:,:,:,times<0),4),[1,1,1,100]) ; 
    bersplow2 = ersplow2 - repmat(mean(ersplow2(:,:,:,times<0),4),[1,1,1,100]) ; 
    bersphigh2 = ersphigh2 - repmat(mean(ersphigh2(:,:,:,times<0),4),[1,1,1,100]) ; 
    
    meanlow1(s,:,:,:) = squeeze(mean(bersplow1(1:2,:,:,:),1)) ; 
    meanhigh1(s,:,:,:) = squeeze(mean(bersphigh1(1:2,:,:,:),1)) ; 
    meanlow2(s,:,:,:) = squeeze(mean(bersplow2(1:2,:,:,:),1)) ; 
    meanhigh2(s,:,:,:) = squeeze(mean(bersphigh2(1:2,:,:,:),1)) ; 
    
    figure,
    subplot(2,2,1) ; imagesc(squeeze(mean(mean(bersplow1,1),2)),[-10,10]) ;
    subplot(2,2,2) ; imagesc(squeeze(mean(mean(bersphigh1,1),2)),[-10,10]) ; 
    subplot(2,2,3) ; imagesc(squeeze(mean(mean(bersplow2,1),2)),[-10,10]) ; 
    subplot(2,2,4) ; imagesc(squeeze(mean(mean(bersphigh2,1),2)),[-10,10]) ; 
end



for s=1:length(subs)  ;
    cd(['c:/shared/badger_eeg/',subs{s},'/outside']) ;  
    low = pop_loadset('epochlow.set') ; 
    high = pop_loadset('epochhigh.set') ; 
    med = pop_loadset('epochmed.set') ; 
    clear ersplow ersphigh
    % for the low frequency filtered data
    for i=1:length(lowcompso{s}) 
        for j=1:15
           [ersplow(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(low.icaact(lowcompso{s}(i),:,j)),low.pnts,[low.xmin,low.xmax],low.srate,0,...
               'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',NaN,'timesout',100,'winsize',64) ;   
        end
    end
    % for the high frequency filtered data
    for i=1:length(highcompso{s}) ;
        for j=1:15
           [ersphigh(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(high.icaact(highcompso{s}(i),:,j)),high.pnts,[high.xmin,high.xmax],high.srate,0,...
               'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',NaN,'timesout',100,'winsize',64) ;  
        end
    end
    otimes = times ; 
    bersplow = ersplow - repmat(mean(ersplow(:,:,:,times<0),4),[1,1,1,100]) ; 
    bersphigh = ersphigh - repmat(mean(ersphigh(:,:,:,times<0),4),[1,1,1,100]) ; 
    meanlow(s,:,:,:) = squeeze(mean(bersplow(1:2,:,:,:),1)) ; 
    meanhigh(s,:,:,:) = squeeze(mean(bersphigh(1:2,:,:,:),1)) ; 
    figure,
    subplot(2,2,1) ; imagesc(squeeze(mean(mean(bersplow,1),2)),[-10,10]) ; 
    subplot(2,2,2) ; imagesc(squeeze(mean(mean(bersphigh,1),2)),[-10,10]) ;
    figure ; 
    for i=1:length(lowcomps{s}) ; subplottight(3,4,i) ; topoplot(squeeze(low.icawinv(:,lowcomps{s}(i))),low.chanlocs) ; colorbar ; end
    for i=1:length(highcomps{s}) ; subplottight(3,4,i+4) ; topoplot(squeeze(high.icawinv(:,highcomps{s}(i))),low.chanlocs) ; colorbar ; end
end


for i=1:8 ; subplot(2,8,i) ;    
    imagesc(freqs,1:15,squeeze(mean(meanhigh1(i,:,:,times>0 & times<3),4)),[-8,8]) ; if i==1 ; xlabel('frequency(hz)') ; ylabel('trial#') ; end
    subplot(2,8,i+8) ; 
    imagesc(freqs,1:15,squeeze(mean(meanlow1(i,:,:,times>0 & times<3),4)),[-8,8]) ; 
    
   % plot(squeeze(mean(meanlow1(i,:,freqs>30 & freqs<100,times>0 & times<5),4))','r') ; hold on ; 
   % plot(squeeze(mean(meanlow2(i,:,freqs>30 & freqs<100,times>0 & times<5),4))','k') ; ylim([-1,9])
end

figure,
for i=1:8 ; subplot(2,8,i) ;    
    %imagesc(freqs,1:15,squeeze(mean(meanhigh1(i,:,:,times>0 & times<3),4)) ,[-8,8]) ; if i==1 ; xlabel('frequency(hz)') ; ylabel('trial#') ; end
    %subplot(2,8,i+8) ; 
    %imagesc(freqs,1:15,squeeze(mean(meanhigh2(i,:,:,times>0 & times<3),4)),[-8,8]) ; 
    
    plot(squeeze(mean(meanhigh1(i,:,freqs>30 & freqs<100,times>0 & times<5),4))','r') ; hold on ; 
    plot(squeeze(mean(meanhigh2(i,:,freqs>30 & freqs<100,times>0 & times<5),4))','k') ; ylim([-1,9]) ;
    if i==1 ; xlabel('frequency(hz)') ; ylabel('modulation (db)') ; ffreqs = freqs(freqs>30 & freqs<100) ; end 
    xlim([0,30]) ; set(gca,'XTick',1:5:length(ffreqs),'XTickLabel',round(ffreqs(1:5:end))) ;
end
subplot(2,8,5) ; plot(1:10,'k') ; hold on ; plot(1:10,'r') ; legend({'stim type 1 (high synchrony)','stim type 2 (low synchrony)'}) ; 

mthigh1 = squeeze(mean(meanhigh1(:,:,:,times>0 & times<5),4)) ; 
mthigh2 = squeeze(mean(meanhigh2(:,:,:,times>0 & times<5),4)) ; 
mtlow1 = squeeze(mean(meanlow1(:,:,:,times>0 & times<5),4)) ; 
mtlow2 = squeeze(mean(meanlow2(:,:,:,times>0 & times<5),4)) ; 

clear tstats ; 
for i=1:8
    [h,p,ci,stats] = ttest2(squeeze(mean(mthigh1(i,:,freqs>40 & freqs<60),3)),squeeze(mean(mthigh2(i,:,freqs>40 & freqs<60),3))) ; 
    tstats(i,1) = stats.tstat ; 
    [h,p,ci,stats] = ttest2(squeeze(mean(mtlow1(i,:,freqs>40 & freqs<60),3)),squeeze(mean(mtlow2(i,:,freqs>40 & freqs<60),3))) ; 
    tstats(i,2) = stats.tstat ; 
end
subplot(1,2,1) ; 
bar(tstats) ; legend({'narrow band ICA','broad band ICA'}) ; xlabel('subject') ; ylabel('t-value') ; 
subplot(1,2,2) ; barwitherr(std(tstats,0,1)./sqrt(8),mean(tstats,1)) ; [h,p,ci,stats] = ttest(tstats(:,1),tstats(:,2)) ; title(['p=',num2str(p)]) ; set(gca,'XTickLabel',{'narrow ICA','broad ICA'}) ; 




errorbar(squeeze(mean(mean(mean(meanhigh(:,:,:,otimes>0 & otimes<3),1),2),4)),squeeze(std(mean(mean(meanhigh(:,:,:,otimes>0 & otimes<3),2),4),0,1))./sqrt(8),'k') ; hold on ; 
errorbar(squeeze(mean(mean(mean(meanhigh1(:,:,:,otimes>0 & otimes<3),1),2),4)),squeeze(std(mean(mean(meanhigh1(:,:,:,otimes>0 & otimes<3),2),4),0,1))./sqrt(8),'m') ; hline(0,'k') ; xlim([0,60]) ; 
set(gca,'XTick',1:5:length(freqs),'XTickLabel',round(freqs(1:5:end))) ;  xlim([0,60]) ; xlabel('frequency (hz)') ; ylabel('modulation (db)') ; legend({'outside','inside'});

for i=1:60 ; [h,p,ci,stats] = ttest(squeeze(mean(mthigh1(:,:,i),2)),squeeze(mean(mean(meanhigh(:,:,i,otimes>0 & otimes<3),2),4))) ; tvals(i) = stats.tstat ;pvals(i) = p ; end
bar(tvals) ;  set(gca,'XTick',1:5:length(freqs),'XTickLabel',round(freqs(1:5:end))) ;  xlim([0,60]) ; xlabel('frequency (hz)') ; ylabel('t-value') ; 
for i=1:length(pvals) ; if pvals(i) < .01 ; text(i,6,'*') ;end ; end ; 
title('* p<0.01') ; 


for i=1:8 ; subplot(2,8,i) ; imagesc(otimes,freqs,squeeze(mean(meanhigh(i,:,:,:),2)),[-8,8]) ; vline([0,3],'k') ; axis xy ; if i==1 ; xlabel('time(s)') ; ylabel('frequency(hz)') ; end ; end
for i=1:8 ; subplot(2,8,i+8) ; imagesc(itimes,freqs,squeeze(mean(meanhigh1(i,:,:,:),2)),[-8,8]) ; axis xy ; vline([0,5],'k') ; end
figure,imagesc([-8,8]) ; colorbar ; 

subplot(1,2,1) ; 
frange = find(freqs>40 & freqs<80) ; 
errorbarxy(squeeze(mean(mean(mean(meanhigh1(:,:,frange,itimes>0 & itimes<3),3),4),2)),squeeze(mean(mean(mean(meanhigh(:,:,frange,itimes>0 & itimes<3),3),4),2)),...
           squeeze(std(mean(mean(meanhigh1(:,:,frange,itimes>0 & itimes<3),3),4),0,2))./sqrt(15),squeeze(std(mean(mean(meanhigh(:,:,frange,itimes>0 & itimes<3),3),4),0,2))./sqrt(15),{'.','k','k'}) ; lsline
[c,p] = corr([squeeze(mean(mean(mean(meanhigh1(:,:,frange,itimes>0 & itimes<3),3),4),2)),squeeze(mean(mean(mean(meanhigh(:,:,frange,itimes>0 & itimes<3),3),4),2))]) ; 
title(['rho=',num2str(c(1,2)),', p=',num2str(p(1,2))]) ; xlabel('gamma(40-80Hz) INSIDE') ; ylabel('gamma(40-80Hz) OUTSIDE') ; 
       
subplot(1,2,2) ; 
frange = find(freqs>7.5 & freqs<25) ; 
errorbarxy(squeeze(mean(mean(mean(meanhigh1(:,:,frange,itimes>0 & itimes<3),3),4),2)),squeeze(mean(mean(mean(meanhigh(:,:,frange,itimes>0 & itimes<3),3),4),2)),...
           squeeze(std(mean(mean(meanhigh1(:,:,frange,itimes>0 & itimes<3),3),4),0,2))./sqrt(15),squeeze(std(mean(mean(meanhigh(:,:,frange,itimes>0 & itimes<3),3),4),0,2))./sqrt(15),{'.','k','k'}) ; lsline
[c,p] = corr([squeeze(mean(mean(mean(meanhigh1(:,:,frange,itimes>0 & itimes<3),3),4),2)),squeeze(mean(mean(mean(meanhigh(:,:,frange,itimes>0 & itimes<3),3),4),2))]) ; 
title(['rho=',num2str(c(1,2)),', p=',num2str(p(1,2))]) ; xlabel('alpha/beta(8-25Hz) INSIDE') ; ylabel('alpha/beta(8-25Hz) OUTSIDE') ; 


for f=1:60
    frange = find(freqs>f*2-2 & freqs<f*2+2) ; 
    [c,p] = corr([squeeze(mean(mean(mean(meanhigh1(:,:,frange,itimes>0 & itimes<3),3),4),2)),squeeze(mean(mean(mean(meanhigh(:,:,frange,itimes>0 & itimes<3),3),4),2))]) ; 
    allc(f) = c(1,2) ; allp(f) = p(1,2) ;   
end



for i=1:8 ; subplot(2,8,i) ; imagesc(itimes,freqs,squeeze(mean(meanlow1(i,:,:,:),2)),[-8,8]) ; vline([0,5],'k') ; axis xy ; if i==1 ; xlabel('time(s)') ; ylabel('frequency(hz)') ; end ; end
for i=1:8 ; subplot(2,8,i+8) ; imagesc(itimes,freqs,squeeze(mean(meanhigh1(i,:,:,:),2)),[-8,8]) ; axis xy ; vline([0,5],'k') ; end
figure,imagesc([-8,8]) ; colorbar ; 


mh1 = squeeze(mean(mean(meanhigh1(:,:,:,itimes>0 & itimes<5),2),4)) ; 
ml1 = squeeze(mean(mean(meanlow1(:,:,:,itimes>0 & itimes<5),2),4)) ; 
clear tstats pstats
for i=1:60 ; [h,p,ci,stats] = ttest(mh1(:,i),ml1(:,i)) ; tstats(i) = stats.tstat ; pstats(i) = p ; end
errorbar(mean(mh1,1),std(mh1,0,1)./sqrt(8),'r') ; xlim([0,60]) ; hline(0,'k') ; hold on ; 
errorbar(mean(ml1,1),std(ml1,0,1)./sqrt(8),'b') ;
set(gca,'XTick',1:5:length(freqs),'XTickLabel',round(freqs(1:5:end))) ;  xlim([0,60]) ; xlabel('frequency (hz)') ; ylabel('modulation(db)') ; legend({'narrow band ICA','broad band ICA'}) ; 
for i=1:length(pstats) ; if pstats(i) < 0.1 ; text(i,4,'*') ; end ; end ; title('* p<0.1') ; 


