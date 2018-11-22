clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie','tegan'} ; 
freqranges = {[48,60],[45,60],[60,80],[50,70],[50,70],[50,70],[56,70],[55,70]} ; 
lowcomps_outside = {[25,6,48],[26,8,18],[15,16,5],[12,5],[9,7],[25,10,3],[13,7],[14,12,11]} ;
highcomps_outside = {[20,17,8],[43,25,50],[48,51],[29,14,24],[63,27],[39,30],[32,25],[51,50,49]} ; 

lowcomps_inside = {[15,25],[35,30,18],[25,17],[35,20,7],[16,15],[14,33],[25,13,35],[37,14,33]} ;
highcomps_inside = {[26,25,21],[44,37],[39,36],[19,38,50],[20,37],[34,19],[54,24,50],[43,49,46]} ; 

% outside
for s=1:length(subs)  ;
    cd(['c:/shared/badger_eeg/',subs{s},'/outside']) ;  
    low = pop_loadset('epochlow.set') ; 
    high = pop_loadset('epochhigh.set') ; 
    clear ersplow ersphigh
    % for the low frequency filtered data
    for i=1:length(lowcomps_outside{s}) 
        for j=1:15
           [ersplow(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(low.icaact(lowcomps_outside{s}(i),:,j)),low.pnts,[low.xmin,low.xmax],low.srate,0,...
               'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',NaN,'timesout',100,'winsize',64) ;   
        end
    end
    % for the high frequency filtered data
    for i=1:length(highcomps_outside{s}) ;
        for j=1:15
           [ersphigh(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(high.icaact(highcomps_outside{s}(i),:,j)),high.pnts,[high.xmin,high.xmax],high.srate,0,...
               'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',NaN,'timesout',100,'winsize',64) ;  
        end
    end
    
    bersplow = ersplow - repmat(mean(ersplow(:,:,:,times<0),4),[1,1,1,100]) ; 
    bersphigh = ersphigh - repmat(mean(ersphigh(:,:,:,times<0),4),[1,1,1,100]) ; 
    
    meanlow_outside(s,:,:,:) = squeeze(mean(bersplow(1:2,:,:,:),1)) ; 
    meanhigh_outside(s,:,:,:) = squeeze(mean(bersphigh(1:2,:,:,:),1)) ; 
    
end
timesout = times ; 

% inside
for s=1:length(subs)  ;
    cd(['c:/shared/badger_eeg/',subs{s}]) ;  
    low = pop_loadset('epochlow.set') ; 
    high = pop_loadset('epochhigh.set') ; 
    clear ersplow ersphigh
    % for the low frequency filtered data
    for i=1:length(lowcomps_inside{s}) 
        for j=1:15
           [ersplow(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(low.icaact(lowcomps_inside{s}(i),:,j)),low.pnts,[low.xmin,low.xmax],low.srate,0,...
               'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',NaN,'timesout',100,'winsize',64) ;   
        end
    end
    % for the high frequency filtered data
    for i=1:length(highcomps_inside{s}) ;
        for j=1:15
           [ersphigh(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(high.icaact(highcomps_inside{s}(i),:,j)),high.pnts,[high.xmin,high.xmax],high.srate,0,...
               'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',NaN,'timesout',100,'winsize',64) ;  
        end
    end
    
    bersplow = ersplow - repmat(mean(ersplow(:,:,:,times<0),4),[1,1,1,100]) ; 
    bersphigh = ersphigh - repmat(mean(ersphigh(:,:,:,times<0),4),[1,1,1,100]) ; 
    
    meanlow_inside(s,:,:,:) = squeeze(mean(bersplow(1:2,:,:,:),1)) ; 
    meanhigh_inside(s,:,:,:) = squeeze(mean(bersphigh(1:2,:,:,:),1)) ; 
    
    impedences = load('impedences.mat') ; 
    imps(s,:,:) = impedences.impedences ; 

end
timesin = times ; 

in = squeeze(mean(meanlow_inside,2)) ; 
out = squeeze(mean(meanlow_outside,2)) ; 

for i=1:60 ; for j=1:100 ; cmat(i,j) = corr2(squeeze(in(:,i,j)),squeeze(out(:,i,j))) ; end ; end

mtin = squeeze(mean(in(:,:,timesin>0 & timesin<5),3)) ; 
mtout = squeeze(mean(out(:,:,timesout>0 & timesout<3),3)) ; 
figure,
plot(squeeze(mean(mtin(:,20:30),2)),squeeze(mean(mtout(:,20:30),2)),'o') ; title(corr2(squeeze(mean(mtin(:,20:30),2)),squeeze(mean(mtout(:,20:30),2)))) ; lsline
xlabel('inside') ; ylabel('outside') ; 
figure,
plot(squeeze(mean(mtin(:,4:8),2)),squeeze(mean(mtout(:,4:8),2)),'o') ; title(corr2(squeeze(mean(mtin(:,4:8),2)),squeeze(mean(mtout(:,4:8),2)))) ; lsline
xlabel('inside') ; ylabel('outside') ; 


meanimps = squeeze(mean(mean(imps,2),3)) ; 
for i=1:60 ; 
    impcorrs(i) = corr2(meanimps,squeeze(mean(mtout(:,i),2))) ; 
    [c,p] = corr([meanimps,squeeze(mean(mtout(:,i),2))]) ;
    subplot(6,10,i) ; plot(meanimps,squeeze(mean(mtout(:,i)-mtin(:,i),2)),'.') ; title(['c=',num2str(c(1,2)),',p=',num2str(p(1,2))]) ; 
end

[c,p] = corr([meanimps,squeeze(mean(mtout(:,4:12),2))]) ;
subplot(1,2,1) ; plot(meanimps,squeeze(mean(mtout(:,4:12),2)),'o') ; lsline ; xlabel('mean impedence') ; ylabel('alpha/beta modulation (8-25Hz)') ; 
title(['rho=',num2str(c(1,2)),',p=',num2str(p(1,2))]) ; 

[c,p] = corr([meanimps,squeeze(mean(mtout(:,15:40),2))]) ;
subplot(1,2,2) ; plot(meanimps,squeeze(mean(mtout(:,15:40),2)),'o') ; lsline ; xlabel('mean impedence') ; ylabel('gamma modulation (30-80Hz)') ; 
title(['rho=',num2str(c(1,2)),',p=',num2str(p(1,2))]) ; 


% plot some spectra (inside vs outside, high vs low) 

errorbar(squeeze(mean(mean(mean(meanlow_outside(:,:,:,timesout>0 & timesout<3),1),2),4)),squeeze(std(mean(mean(meanlow_outside(:,:,:,timesout>0 & timesout<3),2),4),0,1))./sqrt(8)) ; hold on ; 
errorbar(squeeze(mean(mean(mean(meanhigh_outside(:,:,:,timesout>0 & timesout<3),1),2),4)),squeeze(std(mean(mean(meanhigh_outside(:,:,:,timesout>0 & timesout<3),2),4),0,1))./sqrt(8),'r') ; hline(0,'k') ; 
xlim([0,60]) ; legend({'broadband denoising','narrowband denoising'}) ; ylabel('modulation (db)') ; xlabel('frequency(hz)') ; set(gca,'XTick',1:5:length(freqs),'XTickLabel',round(freqs(1:5:end))) ;  

errorbar(squeeze(mean(mean(mean(meanlow_inside(:,:,:,timesin>0 & timesin<3),1),2),4)),squeeze(std(mean(mean(meanlow_inside(:,:,:,timesin>0 & timesin<3),2),4),0,1))./sqrt(8)) ; hold on ; 
errorbar(squeeze(mean(mean(mean(meanhigh_inside(:,:,:,timesin>0 & timesin<3),1),2),4)),squeeze(std(mean(mean(meanhigh_inside(:,:,:,timesin>0 & timesin<3),2),4),0,1))./sqrt(8),'r') ; hline(0,'k') ; 
xlim([0,60]) ; legend({'broadband denoising','narrowband denoising'}) ; ylabel('modulation (db)') ; xlabel('frequency(hz)') ; set(gca,'XTick',1:5:length(freqs),'XTickLabel',round(freqs(1:5:end))) ;  

errorbar(squeeze(mean(mean(mean(meanhigh_outside(:,:,:,timesout>0 & timesout<3),1),2),4)),squeeze(std(mean(mean(meanhigh_outside(:,:,:,timesout>0 & timesout<3),2),4),0,1))./sqrt(8),'k') ; hline(0,'k') ; hold on ; 
errorbar(squeeze(mean(mean(mean(meanhigh_inside(:,:,:,timesin>0 & timesin<3),1),2),4)),squeeze(std(mean(mean(meanhigh_inside(:,:,:,timesin>0 & timesin<3),2),4),0,1))./sqrt(8),'m') ; 
xlim([0,60]) ; legend({'outside narrowband denoised','inside narrowband denoised'}) ; ylabel('modulation (db)') ; xlabel('frequency(hz)') ; set(gca,'XTick',1:5:length(freqs),'XTickLabel',round(freqs(1:5:end))) ;  

% single trial discrimination




