clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','valerie'} ;
comps = {[32,6,48],[9,18,20],[12,10,4],[24,9,16],[8,15,52],[11,21,33],[46,33,15],[17,7,21]} ;
occ = {[65,13],[44,41,28],[60,48],[85,82,67,6],[81,54,36],[77,76,74],[53,49],[77,75,42]} ; 
lingual = {[80,16],[8],[18],[50,5],[60],[53],[31],[16]} ; 
lateral = {[17],[23],[41],[21],[19],[51],[33],[22]} ; 
threshs = [.25,.25,.25,.25,.25,.25,.25,.25] ;
freqs = 1:100 ;  
clear tcorrs ; 
for sub=1:length(subs)
    cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls 
    bcgicas = dir('bcgica*gamma*set') ; 
    
    clear freqepochs boldepochs
    for scan=1:2 
        cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls 
        EEG = pop_loadset(bcgicas(scan).name) ; 
        EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
        allica = load('highfreqs') ; allica = allica.highfreqs ; 
        weights = allica{1} ; sphere = allica{2} ; 
        ica = EEG ; ica.data = weights*sphere*EEG.data ; 
        winv = pinv(weights*sphere) ;
        compacts = ica.data(comps{sub},:) ; 
        fcomps = zeros(size(compacts,1),length(freqs),size(compacts,2)) ; 
        for i=1:length(freqs)
            fcomps(:,i,:) = eegfiltfft(compacts,250,freqs(i)-4,freqs(i)+4) ; 
        end
        mcomps = squeeze(mean(abs(fcomps))) ; 
        events = {ica.urevent.type} ;
        lats = {ica.urevent.latency} ; lats = cell2mat(lats) ; 
        r128s = find(strcmpi('R128',events)) ; 
        startT = round(lats(r128s(1))) ; endT = round(lats(r128s(end))) ;
        stimlats(1,:) = lats((strcmpi('S  1',events)))-startT ; 
        stimlats(2,:) = lats((strcmpi('S  2',events)))-startT ; 
        stimlats(3,:) = lats((strcmpi('S  3',events)))-startT ; 
        scancomps = mcomps(:,startT:endT) ; 
        ntrs = length(r128s) ; 
        rescomps = imresize(scancomps,[size(scancomps,1),ntrs]) ; 
        res_factor = size(scancomps,2)/ntrs ; 
        stimlats = round(stimlats/res_factor) ; 
        
        %fmri
        cd(['c:/shared/newbadger_mri/',subs{sub}]) ;
        corrs = load_untouch_nii('corrs.nii.gz') ; gamma1 = load_untouch_nii(['bp_reg_topup_mc_retino_gamma_0',num2str(scan),'.nii.gz']) ; 
        bincorrs = corrs.img>0.2 ; inds = find(bincorrs==1) ; [cx,cy,cz] = ind2sub(size(corrs.img),inds) ; clear mts ; 
        for i=1:length(inds) ; mts(i,:) = squeeze(gamma1.img(cx(i),cy(i),cz(i),:)) ; end ; mts = mean(mts) ; 
        hrf = spm_hrf(0.693) ; 
        % the resconvs
        padconvs = zeros(size(rescomps,1),size(rescomps,2)+100) ; 
        padconvs(:,101:end) = rescomps ; padconvs(:,1:100) = repmat(mean(rescomps,2),[1,100]) ; 
        for i=1:size(rescomps,1)
            convi = conv(log(padconvs(i,:)),hrf,'full') ;
            padconvs(i,:) = convi(1:size(padconvs,2)) ; 
        end
        padconvs = padconvs(:,101:end) ;
        
        mts(1:10) = 0 ;
        %{
        figure,
        plot(mat2gray(mts),'k','LineWidth',2) ; hold on ; 
        stimcolors = {'b','g','r'} ;
        for i=1:size(stimlats,1)
           for j=1:size(stimlats,2)
              vlines = stimlats(i,j):stimlats(i,j)+7; 
              vline(vlines,stimcolors{i}) ;
           end
        end
        plot(mat2gray(mean(rescomps(40:70,:))),'m','LineWidth',2) ; hold on ; 
        plot(mat2gray(mean(rescomps(10:25,:))),'LineWidth',2,'Color',[0,.7,.7]) ; 
        plot(mat2gray(mts),'k','LineWidth',2) ; xlim([200,500])
        ylabel('signal (A.U.)') ; xlabel('TR') ; 
        plot(1:10,'m','LineWidth',3) ; hold on ; plot(1:10,'LineWidth',3,'Color',[0,.7,.7]) ; plot(1:10,'k','LineWidth',3) ; legend({'40-70hz','10-25Hz','BOLD'})
        %}
        mts = mts(1:size(padconvs,2)) ; 
        for i=1:100
             freqcorrs(sub,scan,i) = corr2(padconvs(i,10:end-10),mts(10:end-10)) ;  
        end
   
        for i=1:size(stimlats,1)
            for j=1:size(stimlats,2)
                freqepochs(scan,i,j,:,:) = padconvs(:,stimlats(i,j):stimlats(i,j)+20) ; 
                substims(sub,scan,i,j,:,:) = rescomps(:,stimlats(i,j)-5:stimlats(i,j)+12) ; 
                boldepochs(scan,i,j,:) = mts(stimlats(i,j):stimlats(i,j)+20) ; 
            end
        end
    end
    allfreqepochs(:,1:16,:,:) = squeeze(freqepochs(1,:,:,:,:)) ; allfreqepochs(:,17:32,:,:) = squeeze(freqepochs(2,:,:,:,:)) ; 
    allboldepochs(:,1:16,:) = squeeze(boldepochs(1,:,:,:)) ; allboldepochs(:,17:32,:) = squeeze(boldepochs(2,:,:,:)) ;     
    subfreqepochs(sub,:,:,:,:) = allfreqepochs ; 

    for i=1:3
        for j=1:32
            for k=1:100
                btcorrs(sub,i,j,k) = corr2(squeeze(allfreqepochs(i,j,k,5:18)),squeeze(allboldepochs(i,j,5:18))) ; 
                
            end
        end
    end
    

    allfreqepochs = allfreqepochs - repmat(mean(allfreqepochs(:,:,:,1),4),[1,1,1,21]) ; 
    mfreqs = squeeze(mean(allfreqepochs(:,:,:,12:16),4)) ;
    mbold = squeeze(mean(allboldepochs(:,:,12:16),3)) ; 
    allboldepochs = allboldepochs - repmat(mean(allboldepochs(:,:,1),3),[1,1,21]) ; 
    
    for i=1:3
        for j=1:100
            for k=1:21
                tcorrs(sub,i,j,k) = corr2(squeeze(mfreqs(i,:,j)),squeeze(allboldepochs(i,:,k))) ; 
            end
        end
    end
    
end


fs = 8:14 ; 
figure,subplot(1,3,1) ; imagesc(squeeze(mean(tcorrs(:,1,:,:),1)),[-.5,.5]);subplot(1,3,2) ; imagesc(squeeze(mean(tcorrs(:,2,:,:),1)),[-.5,.5]);subplot(1,3,3) ; imagesc(squeeze(mean(tcorrs(:,3,:,:),1)),[-.5,.5])
for i=1:100
   [h,p,ci,stats] = ttest(squeeze(mean(tcorrs(:,1,i,fs),4)),squeeze(mean(tcorrs(:,3,i,fs),4))) ; 
   ps(i) = p ; 
   tstats(i) = stats.tstat ; 
   cdiffs(i) = mean((squeeze(mean(tcorrs(:,1,i,fs),4)))-(squeeze(mean(tcorrs(:,3,i,fs),4)))) ; 
end
figure,
subplot(1,2,2) ; 
colors = {'b','g','r'} ;  
shadedErrorBar([],squeeze(mean(mean(tcorrs(:,1,:,fs),1),4)),squeeze(std(mean(tcorrs(:,1,:,fs),4),0,1))/sqrt(8),{'b'}) ; hold on ; hline(0,'k') ;
shadedErrorBar([],squeeze(mean(mean(tcorrs(:,2,:,fs),1),4)),squeeze(std(mean(tcorrs(:,2,:,fs),4),0,1))/sqrt(8),{'g'}) ; hold on ; hline(0,'k') ;
shadedErrorBar([],squeeze(mean(mean(tcorrs(:,3,:,fs),1),4)),squeeze(std(mean(tcorrs(:,3,:,fs),4),0,1))/sqrt(8),{'r'}) ; hold on ; hline(0,'k') ;
ylabel('coupling (r)') ; xlabel('frequency (hz)') ;  
for i=1:length(ps) ; if ps(i) < 0.05 ; text(i,0.35,'*') ; end ; end
title('single trial coupling, *p(0% > 100%) < 0.05, uncorrected') ; 
% legend
subplot(1,2,1) ; 
basestims = (log(substims) - repmat(mean(log(substims(:,:,:,:,:,1:5)),6),[1,1,1,1,1,18]))*10 ; 
meanbase = squeeze(mean(mean(mean(basestims(:,:,:,:,:,7:13),2),4),6)) ; 
shadedErrorBar([],squeeze(mean(meanbase(:,1,:),1)),squeeze(std(meanbase(:,1,:),0,1))/sqrt(8),{'b'}) ; hold on ; 
shadedErrorBar([],squeeze(mean(meanbase(:,2,:),1)),squeeze(std(meanbase(:,2,:),0,1))/sqrt(8),{'g'}) ; 
shadedErrorBar([],squeeze(mean(meanbase(:,3,:),1)),squeeze(std(meanbase(:,3,:),0,1))/sqrt(8),{'r'}) ; hline(0,'k') ; 
xlabel('frequency(hz)') ; ylabel('power(db)') ; title('baseline normalized power') ; 

for i=1:100
   [h,p,ci,stats] = ttest(squeeze(meanbase(:,1,i)),squeeze(meanbase(:,2,i))) ;  
   stimstats(i) = stats.tstat ; 
   stimdiffs(i) = mean(squeeze(meanbase(:,1,i))-squeeze(meanbase(:,3,i))) ; 
end
substimdiffs = squeeze(meanbase(:,1,:) - meanbase(:,3,:)) ; 
subcorrdiffs = squeeze(mean(tcorrs(:,1,:,8:14),4) - mean(tcorrs(:,3,:,8:14),4)) ; 
subplot(2,2,1) ; 
shadedErrorBar([],mean(substimdiffs,1),std(substimdiffs,0,1)/sqrt(8)) ; hline(0,'k') ; xlabel('frequency (hz)') ; ylabel('modulation difference (db)') ; 
subplot(2,2,3) ; 
shadedErrorBar([],mean(subcorrdiffs,1),std(subcorrdiffs,0,1)/sqrt(8)) ;  hline(0,'k') ; xlabel('frequency (hz)') ; ylabel('coupling difference (rho)') ; 
subplot(1,2,2) ; 
plot(stimdiffs,cdiffs,'o') ; lsline
[c,p] = corr(stimdiffs',cdiffs') ; 
title(['r=',num2str(c),', r^2=',num2str(c.^2),', p=',num2str(p)]) ;
xlabel('power difference (db)') ; ylabel('coupling difference (rho)') ; 



mac = squeeze(mean(mean((tcorrs(:,:,:,8:14)),2),4)) ; 
mab = squeeze(mean(meanbase,2)) ; 
subplot(2,2,3) ; shadedErrorBar([],mean(mac,1),std(mac,0,1)./sqrt(8)) ; hline(0,'k') ; xlabel('freq(hz)') ; ylabel('coupling (rho)') ; ylim([-.35,.25])
subplot(2,2,1) ; shadedErrorBar([],mean(mab,1),std(mab,0,1)/sqrt(8)) ; hline(0,'k') ; xlabel('freq(hz)') ; ylabel('power (db)') ; ylim([-3.5,2.5])
subplot(1,2,2) ; plot(mean(mab),mean(mac),'o') ; lsline ; [r,p] = corr(mean(mac)',mean(mab)') ; ylabel('coupling (rho)') ; xlabel('power(db)') ;
title(['r=',num2str(r),', r^2=',num2str(r^2),', p=',num2str(p)]) ; 










