clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 
eegcomps = {[30,26],[32,23],[21,28],[47,46],[19,37],[47,38],[44,46]} ;
epochcomps = {[14,28,30,40],[21,22,45],[15,27,40],[12,15],[13,19],[25,32,36],[33,36,40]} ; 


for sub=1:length(subs)
    %%%% BOLD processing:
    cd(['c:/shared/badger_mri/',subs{sub},'/nii/melodic']) ; 
    mix = load('melodic_mix') ; 
    weights = load_untouch_nii('melodic_IC.nii.gz') ; 
    allweights{sub} = weights.img ; 
    meansub = load_untouch_nii('mean.nii.gz') ; 
    allmeans{sub} = meansub.img ; 
    mixinds = 1:735:735*5 ; 
    for i=1:length(mixinds) ; boldinds{i} = mixinds(i):mixinds(i)+734 ; end
    boldinds{6} = 735*5+1:735*5+1+449 ; 
    for i=1:length(boldinds) ; segmix{i} = mix(boldinds{i},:) ; end
    %bold1 = load_untouch_nii('../reg_topup_mc_retino_gamma_01.nii.gz') ; 
    %bold2 = load_untouch_nii('../reg_topup_mc_retino_gamma_02.nii.gz') ;     
    %boldimgs(1,:,:,:,:) = bold1.img ; boldimgs(2,:,:,:,:) = bold2.img ;     
    %%%% EEG processing:
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    gammas1=dir(['*gamma*Pulse*vhdr']) ; 
    setNames = {gammas1(1).name,gammas1(2).name} ;
    neweeg = pop_loadset('neweeg.set') ; 
    etrigs = {'S  1','S  2','S  3'} ; 
    
    clear trialdata trialTypeLats allEventIndices allbold allersp ; 
    for s=1:length(setNames)
        EEG = pop_loadbv('.',setNames{s}) ; 
        EEG = ica_applyweights(EEG,neweeg) ; 
        %EEG.icaact = eegfiltfft(EEG.icaact,EEG.srate,1,120) ; 
        r128s = find(strcmp('R128',{EEG.urevent.type})) ; 
        lats = {EEG.urevent.latency} ; 
        alltypes = {EEG.urevent.type} ; 
        trlats = cell2mat(lats(r128s)) ; 
        offsetSec = trlats(1)./EEG.srate ; 
        ep = pop_epoch(EEG,{'S  1','S  2','S  3'},[-1,6]) ; 
        ep2 = pop_epoch(EEG,{'S  1','S  2','S  3'},[-7,7]) ; 
        
        types = {ep.epoch.eventtype} ; 
        
        allEventIndices = {ep.epoch.eventurevent} ; 
        eventTypes = {ep.epoch.eventtype} ; 
        trialTypeLats = zeros(3,16) ; 
        trialdata = zeros(3,16,size(ep2.icaact,1),size(ep2.icaact,2)) ;  
        for i=1:length(allEventIndices)
            eventIndicesI = allEventIndices{i} ; 
            for e=1:length(etrigs)
                ind = find(strcmp(etrigs{e},eventTypes{i})) ; 
                if ~isempty(ind) ; break ; end
            end
            lat_i = lats{eventIndicesI{ind}} ;      
            trialdata(e,sum(trialTypeLats(e,:)~=0)+1,:,:) = ep2.icaact(:,:,i) ; 
            trialTypeLats(e,sum(trialTypeLats(e,:)~=0)+1) = lat_i ; 
        end
        
        trialTypeLats = trialTypeLats - trlats(1) ; 
        clear ersp 
        for comp=1:length(epochcomps{sub}) ; 
            for i=1:size(trialdata,1)
                for j=1:size(trialdata,2)
                    [ersp(comp,i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(trialdata(i,j,epochcomps{sub}(comp),:)),ep2.pnts,[ep2.xmin,ep2.xmax],ep2.srate,0,...
                                                                            'plotersp','off','plotitc','off','baseline',NaN,'timesout',100,'freqs',[1,100],'nfreqs',50,'winsize',64) ; 
                end
            end
        end
               
        allersp(:,:,(s-1)*16+1:(s-1)*16+16,:,:) = ersp ; 
        
        smix = segmix{s+2} ; 
        smix = eegfiltfft(smix',1/0.693,0.01,1.5) ; smix = smix' ; 
        trs = round((reshape(trialTypeLats,[1,numel(trialTypeLats)])./EEG.srate)./0.693) ; 
        stimtrs = round(trialTypeLats./(EEG.srate*0.693)) ; 
        design = zeros(1,size(smix,1)) ; stimT = round(5./0.693) ; 
        for i=1:length(trs) ; design(trs(i):trs(i)+stimT) = 1 ; end 
        hrf = spm_hrf(0.693) ; 
        conved = conv(design,hrf,'full') ; 
        conved = conved(1:size(smix,1)) ; 
        
        corrs = corr(smix(50:end-50,:),conved(50:end-50)') ; 
        [sv,si] = sort(corrs,'descend') ; 
        topts = mean(smix(:,si(1:2)),2) ; 
        
        figure,plot(mat2gray(topts(50:end-50))) ; hold on ; plot(mat2gray(conved(50:end-50)),'r') ; legend({'BOLD ICA','glm'}) ; 
        
        baseT = round(1./0.693) ; taskT = round(10./0.693) ; 
        clear boldepochs
        for i=1:size(stimtrs,1)
            for j=1:size(stimtrs,2)
                boldepochs(i,j,:) = topts(stimtrs(i,j)-baseT:stimtrs(i,j)+taskT) - mean(topts(stimtrs(i,j):stimtrs(i,j)+2)) ;    %- squeeze(mean(topts(stimtrs(i,j):topts(stimtrs(i,j)+round(5/0.693))))) ;            
            end
        end  
        allbold(:,(s-1)*16+1:(s-1)*16+16,:) = boldepochs ; 
    end    
   % allersp = allersp - repmat(mean(allersp(:,:,:,:,times<0),5),[1,1,1,1,100 ]) ; 
    mtersp = squeeze(mean(mean(allersp(:,:,:,:,times<5 & times>0),1),5)) ;  
    allmtersp(1:32,:) = squeeze(mtersp(1,:,:)) ; allmtersp(33:64,:) = squeeze(mtersp(2,:,:)) ; allmtersp(65:96,:) = squeeze(mtersp(3,:,:)) ; 
    mbold = squeeze(mean(allbold(:,:,10:end),3)) ; 
    allmbold(1:32) = mbold(1,:) ; allmbold(33:64) = mbold(2,:) ; allmbold(65:96) = mbold(3,:) ; 
    
    for i=1:50 ; allcorrs(sub,i) = corr2(allmtersp(:,i),allmbold') ; end 
    
    for i=1:3 ;
        for j=1:50
            trialcorrs(sub,i,j) = corr2(mbold(i,:),squeeze(mtersp(i,:,j))) ;  
        end
    end 
    allsubersp(sub,:,:,:,:) = squeeze(mean(allersp,1)) ; 
    allsubbold(sub,:,:,:) = allbold ; 
end

figure,
bersp = allsubersp - repmat(mean(allsubersp(:,:,:,:,:,times<0& times>-1),6),[1,1,1,1,1,100]) ; 
meanpow = squeeze(mean(mean(mean(mean(bersp([1,2,3,4,5,6,7],:,:,:,:,times>0 &times<5),1),2),4),6)) ; 
[ax,h1,h2] = plotyy(1:size(meanpow,2),meanpow(1,:),1:size(meanpow,2),mean(allcorrs)) ; hline(0,'k') ; 
set(ax(2),'YTick',-1:.05:1) ; ylabel(ax(2),'single trial correlation (r)') ; ylabel(ax(1),'mean power (db)') ; 
set(ax(1),'XTick',1:5:50,'XTickLabel',round(freqs(1:5:50))) ; set(ax(2),'XTick',1:5:50,'XTickLabel',round(freqs(1:5:50))) ; xlabel('frequency(hz)') ; 
ylim(ax(1),[-2.5,4.5]) ; set(ax(1),'YTick',-2.5:1:4.5) ; 

allbersp = allsubersp - repmat(mean(allsubersp(:,:,:,:,:,times<0& times>-2),6),[1,1,1,1,1,100]) ; 
figure,
for i=1:7 ; subplot(2,7,i) ; 
    plot(squeeze(mean(mean(allbersp(i,:,1,:,15:45,times>0& times<5),2),6))','k') ; hold on ; 
    plot(squeeze(mean(mean(allbersp(i,:,3,:,15:45,times>0& times<5),2),6))','r') ; 
    set(gca,'XTick',1:5:30,'XTickLabel',round(freqs(15:5:45))) ; 
    xlim([1,30]) ; hline(0,'k') ; 
    title(['subject=',num2str(i)]) ; 
    if i==1
    legend({'0%rnd','100%rnd'}) ;  xlabel('frequency(hz)') ; ylabel('power(db)') ; 
    end
end

boldts = -0.693:0.693:10  ; 
for i=1:7 ; subplot(2,7,i+7) ; 
    plot(squeeze(allsubbold(i,1,:,:))','k') ; hold on ; plot(squeeze(allsubbold(i,3,:,:))','r') ; 
    %ylim([min(min(min(squeeze(allsubbold(i,:,:,:))))),max(max(max(squeeze(allsubbold(i,:,:,:)))))]) ;  ; 
    ylim([-6,10]) ;
    if i==1 ;
        
    end
    xlim([1,length(boldts)]) ; set(gca,'XTick',1:3:length(boldts),'XTickLabel',boldts(1:3:end)) ; %vline([1.5/.693,6.5/.693],'k') ;
    hline(0,'k') ; 
    if i==1 ; ylabel('task-baseline') ; xlabel('time(s)') ; end
        title(['subject=',num2str(i)]) ; 

end


mtgamma = squeeze(mean(mean(mean(allbersp(:,:,:,:,freqs>40 & freqs<70,times>0& times<5),2),5),6)) ; 
figure,
z1(:,1:32) = squeeze(mtgamma(:,1,:)) ; z1(:,33:64) = squeeze(mtgamma(:,2,:)) ; z1(:,65:96) = squeeze(mtgamma(:,3,:)) ; 
zref = zeros(1,96) ; zref(1:32) = 1 ; zref(33:64) = 0.9 ;  
for i=1:7 ;
    subplot(2,7,i) ; 
    plot(zref(1:32),z1(i,1:32),'b.'); hold on ;  
    plot(zref(33:64),z1(i,33:64),'r.'); 
    plot(zref(65:96),z1(i,65:96),'k.'); 

    xlim([-.5,1.5]) ; 
    title(['sub=',num2str(i),', r=',num2str(corr2(zref,z1(i,:)))]) ; 
    set(gca,'XTick',0:1,'XTickLabel',{'100%rnd','0%rnd'}) ; 
    if i==1
    ylabel('mean(35-70Hz)') ;
    end
    
    syncorrseeg(i) = corr2(zref,z1(i,:)) ; 
    hline(0,'k') ; 
    
end

zb1(:,1:32) = squeeze(mean(allsubbold(:,1,:,12:15),4)) ; zb1(:,33:64) = squeeze(mean(allsubbold(:,2,:,12:15),4)) ; zb1(:,65:96) = squeeze(mean(allsubbold(:,3,:,12:15),4)) ; 
for i=1:7
    subplot(2,7,i+7) ; 
    plot(zref(1:32),zb1(i,1:32),'b.'); hold on ;  
    plot(zref(33:64),zb1(i,33:64),'r.'); 
    plot(zref(65:96),zb1(i,65:96),'k.'); 
    
    
    xlim([-.5,1.5]) ; 
    title(['sub=',num2str(i),', r=',num2str(corr2(zref,zb1(i,:)))]) ; 
    set(gca,'XTick',0:1,'XTickLabel',{'100%rnd','0%rnd'}) ; 
    if i==1 ; ylabel('BOLD task-rest') ; end
    syncorrs(i) = corr2(zref,zb1(i,:)) ; 
    hline(0,'k') ; 
end
suptitle(['BOLD variability due to synchrony = ',num2str(mean(syncorrs.^2)),'%',', EEG variability due to synchrony = ',num2str(mean(syncorrseeg.^2)),'%']) ; 
%suptitle(['single trial EEG variability explained by synchrony = ',num2str(mean(syncorrseeg.^2)),'%']) ; 

figure,
subplot(1,2,1) ; 
plot(squeeze(mean(mean(mean(allbersp(:,:,1,:,:,times>0 & times<5),2),4),6))','b','LineWidth',2) ; hold on ; 
plot(squeeze(mean(mean(mean(allbersp(:,:,3,:,:,times>0 & times<5),2),4),6))','r','LineWidth',2) ; hline(0,'k') ; 
subplot(1,2,2) ; 
plot(squeeze(mean(allsubbold(:,1,:,:),3))','b','LineWidth',2) ; hold on ; 
plot(squeeze(mean(allsubbold(:,3,:,:),3))','r','LineWidth',2) ; 

figure,
ftimes = times>-6&times<6 ; 
meanfreqs = squeeze(mean(mean(mean(mean(allbersp(:,:,:,:,:,ftimes),2),3),4),6)) ; 
lowf=1 ; highf = 100 ; 
plot(mean(allcorrs(:,freqs>lowf & freqs<highf),2).^2,mean(meanfreqs(:,freqs>lowf & freqs<highf),2).^2,'k.') ; lsline ; 
title(['corr=',num2str(corr(mean(allcorrs(:,freqs>lowf & freqs<highf).^2,2),mean(meanfreqs(:,freqs>lowf & freqs<highf).^2,2)))]) ; 
mtgamma = squeeze(mean(mean(mean(allbersp(:,:,:,:,freqs>35& freqs<70,ftimes),2),5),6)) ; 
z1(:,1:32) = squeeze(mtgamma(:,1,:)) ; z1(:,33:64) = squeeze(mtgamma(:,2,:)) ; z1(:,65:96) = squeeze(mtgamma(:,3,:)) ; 
xlabel('mean r-squared (-6:6s, 1-100Hz)') ; ylabel('modulation 1-100Hz') ; 

% scatter plot the two signals for each subject
figure
for i=1:6 ;
    subplot(2,3,i) ;
    plot(squeeze(zb1(i,1:32)),squeeze(z1(i,1:32)),'b.') ; hold on     
    plot(squeeze(zb1(i,33:64)),squeeze(z1(i,33:64)),'r.') ; 
    plot(squeeze(zb1(i,65:96)),squeeze(z1(i,65:96)),'k.') ; 
    title(num2str(corr2(zb1(i,:),z1(i,:)))) ; 
end

mballsubersp = allsubersp - repmat(mean(allsubersp(:,:,:,:,:,:),6),[1,1,1,1,1,100]) ; 
%mbersp = squeeze(mean(mballsubersp,2)) ; 
mbersp = squeeze(mean(allsubersp,2)) ; 
mbold = squeeze(mean(allsubbold(:,:,:,12:15),4)) ; 
for i=1:7 ;
    for j=1:50
        for k=1:100
            kern(i,j,k) = corr2(squeeze(mbersp(i,:,:,j,k)),squeeze(mbold(i,:,:))) ;    
        end
    end
end

figure,
for i=1:7 ; 
    subplot(1,7,i) ; 
    imagesc(times,freqs,medfilt2(squeeze(kern(i,:,:)),[3,5]),[-.3,.3]) ; vline([0,5],'k') ; axis xy  
    if i==1 ;
       xlabel('time(s)') ; ylabel('freq(hz)') ;   
    end
    title(['sub==',num2str(i)]) ; 
end

figure,
imagesc(times,freqs,medfilt2(squeeze(mean(kern))),[-.2,.2]) ; colorbar ; ylabel('frequency(hz)') ; xlabel('time(s)') ; title('grand avg') ; vline([0,5],'k') ; axis xy ; 


% try to predict the single trials using the kernel. 
mkern = medfilt2(squeeze(mean(kern)),[3,5]) ; 

meeg = squeeze(mean(allsubersp,2)) ; 
for i=1:7
    for j=1:3
        for k=1:32 
            preds(i,j,k) = sum(sum(squeeze(meeg(i,j,k,:,:)).*mkern)) ; 
        end
    end
end

figure
for i=1:7
    subplot(2,4,i) ; hold on ; 
    plot(squeeze(mbold(i,1,:)),squeeze(preds(i,1,:)),'b.') ; stimcorrs(i,1) = corr2(squeeze(mbold(i,1,:)),squeeze(preds(i,1,:))) ; 
    plot(squeeze(mbold(i,2,:)),squeeze(preds(i,2,:)),'r.') ; stimcorrs(i,2) = corr2(squeeze(mbold(i,2,:)),squeeze(preds(i,2,:))) ; 
    plot(squeeze(mbold(i,3,:)),squeeze(preds(i,3,:)),'k.') ; stimcorrs(i,3) = corr2(squeeze(mbold(i,3,:)),squeeze(preds(i,3,:))) ; 
    fullpred(1:32) = squeeze(preds(i,1,:)) ; fullpred(33:64) = squeeze(preds(i,2,:)) ; fullpred(65:96) = squeeze(preds(i,3,:)) ;
    fullbold(1:32) = squeeze(mbold(i,1,:)) ; fullbold(33:64) = squeeze(mbold(i,2,:)) ; fullbold(65:96) = squeeze(mbold(i,3,:)) ; 
    lsline ; 
    if i==1 ; xlabel('BOLD(au)') ; ylabel('EEG(a.u)') ; end
    %r = corr2(fullpred,fullbold) ; 
    [pi,p,ol] = Shepherd(fullpred',fullbold',10000) ; 
    r = corr2(fullpred(ol==0),fullbold(ol==0)) ; 
    title(['subj=',num2str(i),', r=',num2str(r),', r^2=',num2str(r^2)]) ; 
    
end
        
msubraw = squeeze(mean(mean(mean(allsubersp(:,:,1,:,:,times>0 & times<5),2),4),6)) ; 
msubbersp = squeeze(mean(mean(mean(allbersp(:,:,1,:,:,times>0 & times<5),2),4),6)) ; 

msubersp = squeeze(mean(mean(bersp(:,:,:,:,:,:),2),4)) ; 

bersp = allsubersp - repmat(mean(allsubersp(:,:,:,:,:,times>-1 & times<0),6),[1,1,1,1,1,100]) ; 
for i=1:7 ; subplot(3,3,i) ; imagesc(squeeze(mean(mean(bersp(i,:,1,:,:,:),2),4))); end
mbersp = squeeze(mean(mean(mean(bersp(:,:,:,:,:,times>0 & times<5),2),4),6)) ; 


%}





