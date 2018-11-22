clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 
eegcomps = {[30,26],[32,23],[21,28],[47,46],[19,37],[47,38],[44,46]} ;
epochcomps = {[31,32],[13,36],[13,26],[20,40],[14,18],[33,16],[33,37,40]} ; 
compinds = {[16,25,14,12,46,64,80],[8,23,2,31,40,41],[18,40,9,15,48,61],[4,21,1,81,84,67],[53,50,14,74,79,80],[17,20,2,41,75,77],[31,33,31,48,50]} ; % fmri


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
    cd .. ;
    bolds=dir('*bp*gamma*') ; 
    for i=1:length(bolds) ; boldstructs(i) = load_untouch_nii(bolds(i).name) ; end
   
    %boldimgs(1,:,:,:,:) = bold1.img ; boldimgs(2,:,:,:,:) = bold2.img ;     
    %%%% EEG processing:
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    gammas1=dir(['*gamma*Pulse*set']) ; 
    setNames = {gammas1(1).name,gammas1(2).name} ;
    neweeg = pop_loadset('neweeg.set') ; 
    etrigs = {'S  1','S  2','S  3'} ; 
    
    clear trialdata trialTypeLats allEventIndices allbold allersp ; 
    for s=1:length(setNames)
        EEG = pop_loadset(setNames{s}) ; 
        EEG = ica_applyweights(EEG,neweeg) ; 
        %EEG.icaact = eegfiltfft(EEG.icaact,EEG.srate,1,120) ; 
        r128s = find(strcmp('R128',{EEG.urevent.type})) ; 
        lats = {EEG.urevent.latency} ; 
        alltypes = {EEG.urevent.type} ; 
        trlats = cell2mat(lats(r128s)) ; 
        offsetSec = trlats(1)./EEG.srate ; 
        ep = pop_epoch(EEG,{'S  1','S  2','S  3'},[-1,6]) ; 
        ep2 = pop_epoch(EEG,{'S  1','S  2','S  3'},[-7,7]) ; 
        
        % get the continuous time course for visualization purposes
        start = trlats(1)+EEG.srate ; stop = trlats(length(trlats))-EEG.srate ; 
        pows = EEG.icaact(epochcomps{sub},start:stop) ; 
        filtpow = eegfiltfft(pows,EEG.srate,8,15).^2 ; 
        meanfiltpow = squeeze(mean(filtpow,1)) ; 
        resfiltalpha = imresize(meanfiltpow,[1,735]) ; 
        
        filtpow = eegfiltfft(pows,EEG.srate,40,70).^2 ; 
        meanfiltpow = squeeze(mean(filtpow,1)) ; 
        resfiltgamma = imresize(meanfiltpow,[1,735]) ; 
        
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
        smix = eegfiltfft(smix',1/0.693,0.02,1.5) ; smix = smix' ; 
        trs = round((reshape(trialTypeLats,[1,numel(trialTypeLats)])./EEG.srate)./0.693) ; 
        stimtrs = round(trialTypeLats./(EEG.srate*0.693)) ; 
        design = zeros(1,size(smix,1)) ; stimT = round(5./0.693) ; 
        for i=1:length(trs) ; design(trs(i):trs(i)+stimT) = 1 ; end 
        hrf = spm_hrf(0.693) ; 
        conved = conv(design,hrf,'full') ; 
        conved = conved(1:size(smix,1)) ; 
        
        voxcorrs = voxcorr(boldstructs(s).img(:,:,:,40:end-40),conved(40:end-40)) ; 
        goodinds = find(voxcorrs>.4) ;
        [cx,cy,cz] = ind2sub(size(voxcorrs),goodinds) ;
        ts = zeros(length(cx),size(boldstructs(s).img,4)) ; 
        for i=1:length(cx) ; ts(i,:) = squeeze(boldstructs(s).img(cx(i),cy(i),cz(i),:))' ; end
        
        voxcorrs(isnan(voxcorrs)) = 0 ; 
        plotoverlayIntensity2D(meansub.img(:,:,8),mat2gray(abs(voxcorrs(:,:,8))),voxcorrs(:,:,8),270) ; 
        
        corrs = corr(smix(50:end-50,:),conved(50:end-50)') ; 
        [sv,si] = sort(corrs,'descend') ; 
        %topts = mean(smix(:,si(1)),2) ; 
        %topts = mean(smix(:,compinds{sub}(4:5)),2) ; 
        %filttopts = eegfiltfft(topts',1./0.693,0.02,1) ; 
        filttopts = mean(ts,1) ; 
        topts = mean(ts,1) ; 
        
        stimvec = zeros(1,735) ; 
        for i=1:length(trs)
            stimvec(trs(i):trs(i)+round(5./0.693)) = 1 ; 
        end
        
        figure,
        subplot(2,1,1) ; 
        %plot(mat2gray(stimvec(50:250))./2,'k','LineWidth',2) ;
        vline(find(stimvec(50:250)==1),'g') ;
        hold on ;plot(mat2gray(resfiltalpha(50:250)),'LineWidth',2) ;  plot(mat2gray(filttopts(50:250)),'r','LineWidth',2) ;
        plot(mat2gray(resfiltgamma(50:250)),'k','LineWidth',2) ; 
        legend({'alpha(8-15Hz)','BOLD','gamma(40-70Hz)'}) ; ylabel('arb. units') ; xlabel('(TR=0.693)') ; 

        
        
        baseT = round(1./0.693) ; taskT = round(10./0.693) ; 
        clear boldepochs
        for i=1:size(stimtrs,1)
            for j=1:size(stimtrs,2)
                boldepochs(i,j,:) = topts(stimtrs(i,j)-baseT:stimtrs(i,j)+taskT) ;% - mean(topts(stimtrs(i,j)-8:stimtrs(i,j)+2)) ;    %- squeeze(mean(topts(stimtrs(i,j):topts(stimtrs(i,j)+round(5/0.693))))) ;            
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
    
    allsubersp(sub,:,:,:,:) = squeeze(mean(allersp,1)) ; 
    allsubbold(sub,:,:,:) = allbold ; 
end
allsubersp(isnan(allsubersp)) = -40 ; allsubersp(isinf(allsubersp)) = -40 ; 


allsubbersp = allsubersp - repmat(mean(allsubersp(:,:,:,:,times>-8& times<-2),5),[1,1,1,1,100]) ;
mttersp = squeeze(mean(mean(allsubbersp(:,:,:,:,times>0 & times<5),3),5)) ; 
figure,shadedErrorBar([],squeeze(mean(mttersp(:,1,:),1)),squeeze(std(mttersp(:,1,:),0,1))./sqrt(7),'b') ; hold on ; 
shadedErrorBar([],squeeze(mean(mttersp(:,2,:),1)),squeeze(std(mttersp(:,2,:),0,1))./sqrt(7),'g') ; 
shadedErrorBar([],squeeze(mean(mttersp(:,3,:),1)),squeeze(std(mttersp(:,3,:),0,1))./sqrt(7),'r') ; 
hline(0,'k') ; set(gca,'XTick',1:5:50,'XTickLabel',round(freqs(1:5:50))) ; xlabel('frequency(hz)') ;ylabel('power(db)') ; 

boldts = -0.693:0.693:14*0.693 ; 
msubbold = allsubbold - repmat(mean(allsubbold(:,:,:,4:8),4),[1,1,1,16]) ; 
msubbold = squeeze(mean(msubbold,3)) ; 
shadedErrorBar([],squeeze(mean(msubbold(:,1,:),1)),squeeze(std(msubbold(:,1,:),0,1))./sqrt(7),'b') ; hold on ; xlim([1,16]) ; 
shadedErrorBar([],squeeze(mean(msubbold(:,2,:),1)),squeeze(std(msubbold(:,2,:),0,1))./sqrt(7),'g') ; 
shadedErrorBar([],squeeze(mean(msubbold(:,3,:),1)),squeeze(std(msubbold(:,3,:),0,1))./sqrt(7),'r') ; 
set(gca,'XTick',1:2:16,'XTickLabel',boldts(1:2:end)) ; hline(0,'k') ; vline([2,9.2],'r') ; xlabel('time(s)') ; ylabel('BOLD (a.u)') ; 

clear corrs ; 
for i=1:7 ; i 
    for j=1:50
        for k=1:100
            eeg = squeeze(allsubbersp(i,:,:,j,k)) ; eeg2(1:32) = eeg(1,:) ; eeg2(33:64) = eeg(2,:) ; eeg2(65:96) = eeg(3,:) ; 
            bold = squeeze(mean(allsubbold(i,:,:,13:15),4))-squeeze(mean(allsubbold(i,:,:,1:7),4)) ;
            bold2(1:32) = bold(1,:) ; bold2(33:64) = bold(2,:) ; bold2(65:96) = bold(3,:) ; 
            %[r,p] = corr([eeg2;bold2]') ; 
            %corrs(i,j,k) = r(1,2) ; 
            %ps(i,j,k) = p(1,2) ; 
            mdl = LinearModel.fit(eeg2,bold2) ; 
            corrs(i,j,k) = mdl.Coefficients.tStat(2) ; 
            pvals(i,j,k) = mdl.Coefficients.pValue(2) ; 
        end
    end
end


clear fcorrs ; 
for i=1:7 ; i 
    for j=1:50
            eeg = squeeze(mean(allsubbersp(i,:,:,j,times>0 & times<5),5)) ; eeg2(1:32) = eeg(1,:) ; eeg2(33:64) = eeg(2,:) ; eeg2(65:96) = eeg(3,:) ; 
            bold = squeeze(mean(allsubbold(i,:,:,13:15),4))-squeeze(mean(allsubbold(i,:,:,1:7),4)) ;
            bold2(1:32) = bold(1,:) ; bold2(33:64) = bold(2,:) ; bold2(65:96) = bold(3,:) ; 
            %[r,p] = corr([eeg2;bold2]') ; 
            %corrs(i,j,k) = r(1,2) ; 
            %ps(i,j,k) = p(1,2) ; 
            mdl = LinearModel.fit(eeg2,bold2) ; 
            fcorrs(i,j) = mdl.Coefficients.tStat(2) ; 
            fpvals(i,j) = mdl.Coefficients.pValue(2) ; 
    end
end


plot(fcorrs','k') ; hline(0,'k') ; hold on ; shadedErrorBar([],mean(fcorrs,1),std(fcorrs,0,1)./sqrt(7),'r') ;
set(gca,'XTick',1:5:50,'XTickLabel',round(freqs(1:5:50))) ; xlabel('frequency(hz)') ;ylabel('t-value') ; 
for i=1:50
   [h,p,ci,stats] = ttest(fcorrs(:,i)) ;  
    ps(i) = p ; 
    if ps(i) < 0.01 ; text(i,5.5,'*') ; end
end
title('*p<0.01, uncorrected') ; 




mtcorrs = squeeze(mean(corrs(:,:,times>0 & times<5),3)) ; 
mttersp = squeeze(mean(mean(mean(allsubbersp(:,:,:,:,times>0 & times<5),2),3),5)) ; 

plot(mean(mtcorrs(:,20:35),2),mean(mttersp(:,20:35),2),'o') ; lsline ; 
[c,p] = corr([mean(mtcorrs(:,20:35),2),mean(mttersp(:,20:35),2)]) ; title(['rho=',num2str(c(1,2)),', p=',num2str(p(1,2))]) ; 
xlabel('40-70Hz coupling (t-value)') ; ylabel('power modulation (db)') ; xlim([-.2,1.8]) ; 


figure, 
plot(squeeze(mean(mttersp,1)),squeeze(mean(mtcorrs,1)),'o') ; lsline ; 
[c,p] = corr([squeeze(mean(mttersp,1))',squeeze(mean(mtcorrs,1))']) ; 
title(['rho=',num2str(c(1,2)),', p=',num2str(p(1,2))]) ; xlim([-1.5,1.5]) ; xlabel('power(db)') ; ylabel('coupling(t)') ; 
subplot(1,2,2) ; 
shadedErrorBar([],mean(mttersp,1),std(mttersp,0,1)./sqrt(7),'b') ; hold on ; 
shadedErrorBar([],mean(mtcorrs,1),std(mtcorrs,0,1)./sqrt(7),'r') ; 


imagesc(times,freqs,medfilt2(squeeze(mean(corrs,1)))) ; axis xy ;
coupling = medfilt2(squeeze(mean(corrs,1))) ; 
plot(mean(coupling(:,times>0 & times<5),2)) ; hline(0,'k') ; 
titles = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','GRAND AVG'} ; 
for i=1:8  ;  
    if i<8
        subplot(1,8,i) ; imagesc(times,freqs,medfilt2(squeeze(corrs(i,:,:))),[-4,4]) ; 
        axis xy ; vline([0,5],'k') ;
        title(titles(i)) ; 
        meds(i,:,:) = medfilt2(squeeze(corrs(i,:,:))) ; 
    else
        subplot(1,8,i) ; imagesc(times,freqs,squeeze(mean(meds,1)),[-2,2]) ; 
        axis xy ; vline([0,5],'k') ;
        title(titles(i)) ; 
    end
end 

%shadedErrorBar([],squeeze(mean(mean(corrs(:,:,times>0 & times<5),1),3)),squeeze(std(mean(corrs(:,:,times>0 & times<5),3),0,1))./sqrt() ; 










allsubbersp = allsubersp - repmat(mean(allsubersp(:,:,:,:,times>-7& times<-2),5),[1,1,1,1,100]) ;
colors = {'b','g','r'} ;
clear meancorrs ; 
clear stcorrs stpvals
for x=1:3
for i=1:7 ; 
    for j=1:50
        meancorrs(x,i,j) = corr2(squeeze(mean(allsubbersp(i,x,:,j,times>0 & times<5),5)),squeeze(mean(allsubbold(i,x,:,13:15),4))-squeeze(mean(allsubbold(i,x,:,5:7),4))) ; 
        
        mdl = LinearModel.fit(squeeze(mean(allsubbersp(i,x,:,j,times>0 & times<5),5)),squeeze(mean(allsubbold(i,x,:,13:15),4))-squeeze(mean(allsubbold(i,x,:,5:7),4))) ; 
        stcorrs(x,i,j) = mdl.Coefficients.tStat(2) ; 
        stpvals(x,i,j) = mdl.Coefficients.pValue(2) ; 
        
    end
end
end


for j=1:50
    [h,p1(j),ci,stats] = ttest(squeeze(stcorrs(1,:,j)),squeeze(stcorrs(2,:,j))) ; 
    [h,p2(j),ci,stats] = ttest(squeeze(stcorrs(1,:,j)),squeeze(stcorrs(3,:,j))) ; 
    [h,p3(j),ci,stats] = ttest(squeeze(stcorrs(2,:,j)),squeeze(stcorrs(3,:,j))) ; 
    [h,p4(j),ci,stats] = ttest(squeeze(mean(stcorrs(1:2,:,j),1)),squeeze(stcorrs(3,:,j))) ; 
    
end

shadedErrorBar([],squeeze(mean(mean(stcorrs(1:2,:,:),1),2)),squeeze(std(mean(stcorrs(1:2,:,:),1),0,2))./sqrt(7),'b') ; hold on ; 
shadedErrorBar([],squeeze(mean(mean(stcorrs(3,:,:),1),2)),squeeze(std(mean(stcorrs(3,:,:),1),0,2))./sqrt(7),'r') ;
ylim([-2.5,2.5]) ; hline(0,'k') ; 
for i=1:50
    if p4(i) < 0.05 ; text(i,2.2,'*') ; end
end
set(gca,'XTick',1:5:50,'XTickLabel',round(freqs(1:5:50))) ; xlabel('frequency(hz)') ;ylabel('t-value') ; title('*p<0.05,uncorrected') ; 









for i=1:50 ; 
   [h,p,ci,stats] = ttest(squeeze(mean(meancorrs(1:2,:,i),1)),squeeze(meancorrs(3,:,i))) ; 
   allts(i) = stats.tstat ; 
   allps(i) = p ; 
end
%shadedErrorBar([],squeeze(mean(meancorrs([1],:,:),2))',squeeze(std(meancorrs([1],:,:),0,2))'./sqrt(7),'b') ; 
shadedErrorBar([],squeeze(mean(mean(meancorrs([1:2],:,:),2),1))',squeeze(std(mean(meancorrs([1:2],:,:),1),0,2))'./sqrt(7),'b') ;hold on ; 
shadedErrorBar([],squeeze(mean(meancorrs([3],:,:),2))',squeeze(std(meancorrs([3],:,:),0,2))'./sqrt(7),'r') ; 
hline(0,'k') ; xlabel('frequency(hz)') ; ylabel('correlation (r)') ; hold on ; 
set(gca,'XTick',1:5:50,'XTickLabel',1:10:100) ; xlim([1,50])
for i=1:length(allps) ; if allps(i) < 0.05 ; text(i,0.3,'*');  end ; end

mtersp = squeeze(mean(allsubbersp(:,:,:,:,times>-1 & times<5),5)) ; 
mbold = squeeze(mean(allsubbold(:,:,:,13:15),4))-squeeze(mean(allsubbold(:,:,:,5:7),4)) ; 
for i=1:7 
   subplot(2,8,i) ; 
   boldi = reshape(squeeze(mbold(i,:,:)),[1,96]) ; 
   eegi = reshape(squeeze(mean(mtersp(i,:,:,20:30),4)),[1,96]) ; %- reshape(squeeze(mean(mtersp(i,:,:,1:15),4)),[1,96]) ; 
   plot(boldi,eegi,'m.') ; lsline ; [c,p] = corr([boldi;eegi]') ; 
   title(['rho=',num2str(c(1,2)),' p=',num2str(p(1,2))]) ; if i==1 ; xlabel('bold') ; ylabel('gamma') ; end
   
   subplot(2,8,i+8) ; 
   boldi = reshape(squeeze(mbold(i,:,:)),[1,96]) ; 
   eegi = reshape(squeeze(mean(mtersp(i,:,:,4:7),4)),[1,96]) ; %- reshape(squeeze(mean(mtersp(i,:,:,1:15),4)),[1,96]) ; 
   plot(boldi,eegi,'.') ; lsline ; [c,p] = corr([boldi;eegi]') ; 
   title(['rho=',num2str(c(1,2)),' p=',num2str(p(1,2))]) ; if i==1 ; xlabel('bold') ; ylabel('alpha') ; end
end



[c,p] = corr([squeeze(mean(mean(mean(mtersp(:,:,:,20:30),2),3),4)),squeeze(mean(meancorrs(:,20:30),2))]) ; 
plot(squeeze(mean(mean(mean(mtersp(:,:,:,20:30),2),3),4)),squeeze(mean(meancorrs(:,20:30),2)),'o') ; lsline
title(['rho=',num2str(c(1,2)),' p=',num2str(p(1,2))]) ; xlabel('mean gamma power') ; ylabel('mean gamma coupling (r)') ; 

subplot(1,2,1) ; imagesc([-.35,.35]) ; colorbar ;  subplot(1,2,2) ; imagesc([-.2,.2]) ; colorbar ;  


errorbar(squeeze(mean(mean(mean(allsubbersp(:,:,:,:,times>0 & times<5),1),3),5))',squeeze(std(mean(mean(allsubbersp(:,:,:,:,times>0 & times<5),5),3),0,1))'./sqrt(7)) ; hline(0,'k') ; xlim([0,51]) ;
set(gca,'XTick',1:5:length(freqs),'XTickLabel',round(freqs(1:5:end))) ; xlabel('frequency(hz)') ; ylabel('power(db)') ; legend({'0%rnd','10%rnd','100%rnd'}) ; 

basenormbold = allsubbold - repmat(mean(allsubbold(:,:,:,5:7),4),[1,1,1,16]) ;
errorbar(squeeze(mean(mean(basenormbold,1),3))',squeeze(std(mean(basenormbold,3),0,1))'./sqrt(7)) ; 
stimT = 5./0.693 ; vline([baseT+1,baseT+stimT+1]) ; xt = -0.693:0.693:0.693*14 ; set(gca,'XTick',1:3:16,'XTickLabel',xt(1:3:end)) ; xlabel('time(s)') ; ylabel('BOLD (A.U)') ; 

for i=1:7 ; subplot(1,7,i) ; imagesc(times,freqs,squeeze(mean(mean(allsubbersp(i,:,:,:,:),2),3)),[-5,5]) ; vline([0,5],'k') ;  axis xy ; if i==1 ; xlabel('time(s)') ; ylabel('freq(hz)') ;end ; title(['sub ',num2str(i)]) ; end


for i=1:7 ; 
    subplot(2,4,i) ;
    pre = squeeze(reshape(allsubbold(i,:,:,6),[1,3*32])) ; peak = squeeze(reshape(allsubbold(i,:,:,14),[1,3*32])) ; 
    plot(pre,peak,'.') ;  
    title(['corr2=',num2str(corr2(pre,peak))]) ; lsline ; 
    ylabel('peak (next trial)') ; xlabel('undershoot (prev. trial)') ; 
end
suptitle('BOLD peak as a function of previous trial undershoot') ; 

lowf = 8 ; highf = 15 ; 
for i=1:7
   subplot(2,4,i) 
   prepow = squeeze(reshape(mean(mean(allsubersp(i,:,:,freqs>lowf & freqs<highf,times>-7 & times<-3),4),5),[1,96])) ; 
   postpow = squeeze(reshape(mean(mean(allsubersp(i,:,:,freqs>lowf & freqs<highf,times>0 & times<5),4),5),[1,96])) ; 
   plot(prepow,postpow,'.') ; lsline ; title(['corr2=',num2str(corr2(prepow,postpow))]) ; 
   
end



% power dependence:
subplot(1,2,1) ; 
plot(mat2gray(squeeze(mean(mean(corrs(:,:,times>0 & times<5),1),3))));hold on;plot(mat2gray(squeeze(mean(mean(mean(mean(allsubbersp(:,:,:,:,times>0 & times<5),1),2),3),5))),'r') ; lsline
title(corr2(mat2gray(squeeze(mean(mean(corrs(:,:,times>0 & times<5),1),3)))',mat2gray(squeeze(mean(mean(mean(mean(allsubbersp(:,:,:,:,times>0 & times<5),1),2),3),5))))) ; 

subplot(1,2,2) ; 
plot((squeeze(mean(mean(corrs(:,:,times>0 & times<5),1),3))),(squeeze(mean(mean(mean(mean(allsubbersp(:,:,:,:,times>0 & times<5),1),2),3),5))),'.') ; lsline ; hline(0,'k-') ; vline(0,'k-') ; 
title(corr2((squeeze(mean(mean(corrs(:,:,times>0 & times<5),1),3)))',(squeeze(mean(mean(mean(mean(allsubbersp(:,:,:,:,times>0 & times<5),1),2),3),5))))) ; 


% synchrony analysis:
synch = [0,.1,1] ; 
repsynch = (repmat(synch,[32,1])) ; synch1 = [repsynch(:,1);repsynch(:,2);repsynch(:,3)] ; 
for s=1:7 ; subplot(1,7,s) ; 

for i=1:50 ; subplot(5,10,i) ; 
    plot(ones(1,32)*0,squeeze(mean(mean(allsubersp(s,1,:,freqs>i*2-5 & freqs<i*2+5,times>0 & times<5),4),5)),'b.') ; hold on ; 
    plot(ones(1,32)*.1,squeeze(mean(mean(allsubersp(s,2,:,freqs>i*2-5 & freqs<i*2+5,times>0 & times<5),4),5)),'g.') ; 
    plot(ones(1,32)*1,squeeze(mean(mean(allsubersp(s,3,:,freqs>i*2-5 & freqs<i*2+5,times>0 & times<5),4),5)),'r.') ; xlim([-.1,1.1]) ; 
    a(1:32) = squeeze(mean(mean(allsubersp(s,1,:,freqs>i*2-5 & freqs<i*2+5,times>0 & times<5),4),5)) ; 
    a(33:64) = squeeze(mean(mean(allsubersp(s,2,:,freqs>i*2-5 & freqs<i*2+5,times>0 & times<5),4),5)) ; 
    a(65:96) = squeeze(mean(mean(allsubersp(s,3,:,freqs>i*2-5 & freqs<i*2+5,times>0 & times<5),4),5)) ; 
    allsyncs(s,i) = (corr2(synch1',a)) ; 
end
    
    a(1:32) = squeeze(mean(mean(allsubersp(s,1,:,freqs>40 & freqs<70,times>0 & times<5),4),5)) ; 
    a(33:64) = squeeze(mean(mean(allsubersp(s,2,:,freqs>40 & freqs<70,times>0 & times<5),4),5)) ; 
    a(65:96) = squeeze(mean(mean(allsubersp(s,3,:,freqs>40& freqs<70,times>0 & times<5),4),5)) ; 
    fsynchs(s) = (corr2(synch1',a)) ; 
    plot(synch1,a,'.') ; xlim([-.2,1.2]) ; lsline ; title(['r=',num2str(fsynchs(s))]) ; if s==1 ; ylabel('log power (40-70Hz') ; end
    set(gca,'XTick',[0,.1,1],'XTickLabel',{'0','10','100'}) ; xlabel('%SR') ; 
end


shadedErrorBar([],mean(allsyncs.^2),std(allsyncs,0,1)./sqrt(7)) ; suptitle('variance explained by synchrony') ; 
set(gca,'XTick',1:5:50,'XTickLabel',round(freqs(1:5:50))) ; ylabel('r^2') ; xlabel('frequency(hz)') ; hline(0,'k') ; 

[c,p] = corr([mean(allsyncs)',mean(mttersp)']) ; 
plot(mean(allsyncs),mean(mttersp),'o') ; lsline
xlabel('synchrony dependence (rho)') ; ylabel('power modulation (db)') ; title('rho=',num2str(c(1,2)),', p=',num2str(p(1,2))]) ; 

% BOLD synchrony analysis:
mbold = squeeze(mean(allsubbold(:,:,:,2:7),4) - mean(allsubbold(:,:,:,12:15),4)) ; 
for i=1:7
    a(1:32) = mbold(i,1,:) ; a(33:64) = mbold(i,2,:) ; a(65:96) = mbold(i,3,:) ; 
    bcorrs(i) = corr2(a',synch1) ; 
    subplot(1,7,i) ; plot(synch1,a,'.') ; lsline ; title(['r=',num2str(bcorrs(i))]) ; xlim([-.2,1.2]) ; set(gca,'XTick',[0,.1,1],'XTickLabel',{'0','10','100'}) ; xlabel('%SR') ; if i==1 ; ylabel('BOLD peak') ; end
end

barwitherr(std(bcorrs.^2)./sqrt(7),mean(bcorrs.^2)) ; ylim([-.1,1]) ; xlim([.5,1.5]) ; set(gca,'XTick',[]) ; xlabel('BOLD') ; ylabel('r^2') ; 






