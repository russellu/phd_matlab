clear all ; close all ; 
subs = {'lowcont_russell'} ; 
comps = {[7,13,15]} ; 
stims = {'S  2','S  4'} ;

for sub=1:length(subs)
    
    cd(['E:/jly_orientation/',subs{sub},'/or']) ;  
    discr = dir('*vhdr') ; 
    for i=1:length(discr)
       EEG = pop_loadbv('.',discr(i).name) ; 
       EEG = pop_resample(EEG,256) ; 
       if i==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG) ; end
    end
    merged = pop_chanedit(merged,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
    merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61) - eegfiltfft(merged.data,merged.srate,84,86) ;  
    
    %figure,bar(sum(abs(diff(merged.data,1,2)),2))
    merged = pop_interp(merged,[16,62],'spherical');
    
    %{
    
    mergefilt = merged ; mergefilt.data = eegfiltfft(merged.data,merged.srate,1,90)  ; 
    [weights,sphere] = runica(mergefilt.data(:,1:3:end),'maxsteps',128); 
    ica = merged; 
    ica.icaact = weights*sphere*merged.data ; winv = pinv(weights*sphere) ; 
    icaw{1} = weights; icaw{2} = sphere ; save('icaw','icaw') ; 
    ica.data = ica.icaact ; clear allersp ; 
    newmerged = merged ; newmerged.data = weights*sphere*merged.data ; 
    for s=1:length(stims) ; disp(s) ; 
        allep = pop_epoch(newmerged,{stims{s}},[-1,21]) ; 
    for i=1:64 
            [allersp(s,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,14:end)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',0,'verbose','off','timesout',200) ; 
    end
    end
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(allersp(1,i,:,:)),[-8,8]) ; colormap jet ;  title(i) ; axis xy ; end
    %}
        
    icaw = load('icaw') ; icaw = icaw.icaw ; winv = pinv(icaw{1}*icaw{2}) ; 
    newmerged = merged ; newmerged.data = icaw{1}*icaw{2}*merged.data ; 
    goodcs = comps{sub} ; 
    clear stersp 
    for s=1:length(stims) ; disp(s) ; 
        allep = pop_epoch(newmerged,{stims{s}},[-2,19]) ; 
    for i=1:length(goodcs)
        for k=1:size(allep.data,3)
            [stersp(s,i,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(goodcs(i),:,k)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',NaN,'verbose','off','timesout',200) ; 
        end
    end
    end
    bstersp = stersp - repmat(mean(stersp(:,:,:,:,times<0 | times>18),5),[1,1,1,1,200]) ; 
    mbstersp = squeeze(mean(bstersp,2)) ; 
  
    subts(1,:,:,:) = squeeze(mean(mbstersp(:,1:30,:,:),2)); 
    subts(2,:,:,:) = squeeze(mean(mbstersp(:,31:60,:,:),2)); 
    subts(3,:,:,:) = squeeze(mean(mbstersp(:,61:90,:,:),2)); 

    %{
    subersp{sub} = bstersp ; 
    figure, 
    for i=1:2 
        subplot(2,2,i) ; goodts = zeros(1,size(mbstersp,2)) ; 
        goodts(badts{sub}{i}) = 1 ; goodts = find(goodts==0) ; colormap jet ; axis xy ; 
        imagesc(squeeze(mean(mbstersp(i,goodts,:,:),2)),[-4,4]) ; 
        subts(sub,i,:,:) = squeeze(mean(mbstersp(i,goodts,:,:),2)) ; 
    end
    %}
end


revrots = load('E:\jly_orientation\revrots') ; revrots = revrots.revrots ; 
fwdrots = load('E:\jly_orientation\fwdrots') ; fwdrots = fwdrots.fwdrots ; 
stimes = find(times>=0 & times<=18) ; 
revrots = imresize(revrots,[1,length(stimes)]) ; 
fwdrots = imresize(fwdrots,[1,length(stimes)]) ; 
[sv,rsi] = sort(revrots,'ascend') ; 
[~,fsi] = sort(fwdrots,'ascend') ; 
subplot(1,2,2) ; 
imagesc(sv,freqs,squeeze(mean(subts(:,1,:,stimes(fsi)),1)) + squeeze(mean(subts(:,2,:,stimes(rsi)),1)),[-6,6]) ; axis xy ; colormap jet ; 
xlabel('bar orientation (deg) (vertical bars = 0)') ; ylabel('frequency(hz)') ; 

rets = squeeze((subts(:,1,:,stimes(fsi)))) + squeeze((subts(:,3,:,stimes(rsi)))) ; 
orients = squeeze((subts(:,1,:,stimes(fsi)))) + squeeze((subts(:,2,:,stimes(rsi)))) ; 
for i=1:4 ; retcorrs(i,:,:) = corr(squeeze(rets(i,:,:))') ; end ; 
[meanretcorrs,p] = corr(squeeze(mean(rets,1))') ; [meanorientcorrs,p2] = corr(squeeze(mean(orients,1))') ;
subplot(1,2,1) ; imagesc(freqs,freqs,meanretcorrs,[-1,1]) ; colormap jet ; xlabel('frquency(hz)') ; ylabel('frequency(hz)') ; title('retinotopic tuning correlations') ; 
subplot(1,2,2) ; imagesc(freqs,freqs,meanorientcorrs,[-1,1]) ; colormap jet ; xlabel('frquency(hz)') ; ylabel('frequency(hz)') ; title('orientation tuning correlations') ; 

clear fwdinds revinds
icount = 1 ; 
for i=1:400
   fwdinds{icount} = find(fwdrots<i+10 & fwdrots>i-10) ; 
   revinds{icount} = find(revrots<i+10 & revrots>i-10) ; 
   icount = icount + 1 ; 
end
clear orvals retvals ; 
for i=1:length(fwdinds)
   orvals(:,:,:,i) = mean(subts(:,1,:,stimes(fwdinds{i})),4) + mean(subts(:,2,:,stimes(revinds{i})),4) ; 
end

titles = {'5%','15%','100%'};
for i=1:3 ; subplot(2,3,i) ; imagesc(squeeze(orvals(i,1,:,:)),[-4,4]) ; colormap jet; axis xy ; title(titles{i}); end ; xlabel('orientation(deg)'); 
subplot(2,3,4) ;
plot(squeeze(mean(orvals(2,1,14:16,:),3)));  hold on ; plot(squeeze(mean(orvals(3,1,14:16,:),3))); title('28-32Hz tuning (low gamma)'); xlabel('orientation(deg)'); 
legend({'15% contrast', '100% contrast'});
subplot(2,3,5) ;
plot(squeeze(mean(orvals(2,1,22:30,:),3)));  hold on ; plot(squeeze(mean(orvals(3,1,22:30,:),3))); title('40-60Hz tuning (high gamma)'); xlabel('orientation(deg)'); 
legend({'15% contrast', '100% contrast'});
%{
for i=1:length(subersp)
    subi = squeeze(subersp{i}) ; ntrials = size(subi,3) ; 
    clear suborvals subretvals
    for j=1:length(fwdinds)
        suborvals(1:ntrials,:,j) = squeeze(mean(mean(subi(2,:,:,:,stimes(fwdinds{j})),5),2)) ; 
        suborvals(ntrials+1:ntrials*2,:,j) = squeeze(mean(mean(subi(4,:,:,:,stimes(revinds{j})),5),2)) ;       
        subretvals(1:ntrials,:,j) = squeeze(mean(mean(subi(1,:,:,:,stimes(fwdinds{j})),5),2)) ; 
        subretvals(ntrials+1:ntrials*2,:,j) = squeeze(mean(mean(subi(3,:,:,:,stimes(revinds{j})),5),2)) ;         
    end
    allsuborvals{i} = suborvals ; 
    allsubretvals{i} = subretvals ; 
end
nsubs = length(subs) ; 

for i=1:length(allsuborvals) ; subplot(1,nsubs+1,i) ;
    badis = badts{i} ; ntrials = size(allsuborvals{i},1)/2 ; 
    badis1 = badis{1} ; badis2 = badis{3}+ntrials ; goods = zeros(1,(ntrials)*2) ; goods([badis1,badis2]) = 1 ; goods = find(goods==0) ; 
    imagesc(squeeze(mean(allsubretvals{i}(goods,:,:),1)),[-2,2]) ; axis xy ; colormap hot ; 
    allrets(i,:,:) = squeeze(mean(allsubretvals{i}(goods,:,:),1)) ; 
    orrsi = squeeze((allsubretvals{i}(goods,:,50:end-50))) ; 
    for j=1:60
       meanij = squeeze(mean(orrsi(:,j,:),1)) ; 
       stdij = squeeze(std(orrsi(:,j,:),0,1)) ; 
       maxij = find(meanij==max(meanij)) ; minij = find(meanij==min(meanij)) ; 
       rdprimes(i,j) = (meanij(maxij)-meanij(minij)) ./ (0.5*(stdij(maxij)+stdij(minij))) ;      
    end
end 
subplot(1,nsubs+1,nsubs+1) ; imagesc(squeeze(mean(allrets,1)),[-1,1]) ; colormap hot ; axis xy ; 

% compute dprime
for i=1:length(allsuborvals) ; subplot(1,nsubs,i) ;
    badis = badts{i} ; ntrials = size(allsuborvals{i},1)/2 ; 
    badis1 = badis{2} ; badis2 = badis{4}+ntrials ; goods = zeros(1,(ntrials)*2) ; goods([badis1,badis2]) = 1 ; goods = find(goods==0) ; 
    imagesc(squeeze(mean(allsuborvals{i}(goods,:,:),1)),[-6,6]) ; axis xy ; colormap jet ; 
    orrsi = squeeze((allsuborvals{i}(goods,:,50:end-50))) ; 
    for j=1:60
       meanij = squeeze(mean(orrsi(:,j,:),1)) ; 
       stdij = squeeze(std(orrsi(:,j,:),0,1)) ; 
       maxij = find(meanij==max(meanij)) ; minij = find(meanij==min(meanij)) ; 
       odprimes(i,j) = (meanij(maxij)-meanij(minij)) ./ (0.5*(stdij(maxij)+stdij(minij))) ; 
        
    end
    allors(i,:,:) = squeeze(mean(allsuborvals{i}(goods,:,:),1)) ; 
end 

for i=1:60 
    for j=1:60
        [h,p,ci,stats] = ttest(odprimes(:,i),odprimes(:,j)) ; 
        ostats(i,j) = stats.tstat ; ops(i,j) = p ; 
        [h,p,ci,stats] = ttest(rdprimes(:,i),rdprimes(:,j)) ; 
        rstats(i,j) = stats.tstat ; rps(i,j) = p ; 
    end
end
subplot(2,2,2) ; 
shadedErrorBar([],mean(odprimes),std(odprimes,0,1)./sqrt(5),{'g'}) ; hold on ;
shadedErrorBar([],mean(rdprimes),std(rdprimes,0,1)./sqrt(5),{'k'}) ;  title('grand average d-prime') ; 
set(gca,'XTick',1:5:60,'XTickLabel',[1,(10:10:120)]) ; xlabel('frequency(hz)') ; ylabel('tuning (d-prime)') ; 
subplot(2,2,1) ; 
shadedErrorBar([],squeeze(mean(mean(allors(:,:,50:350),1),3)),squeeze(std(mean(allors(:,:,50:350),3),0,1))./sqrt(5),{'g'}) ; hold on ; 
shadedErrorBar([],squeeze(mean(mean(allrets(:,:,50:350),1),3)),squeeze(std(mean(allrets(:,:,50:350),3),0,1))./sqrt(5),{'k'}) ; hold on ; hline(0,'k') ;
title('grand average power modulation') ; set(gca,'XTick',1:5:60,'XTickLabel',[1,(10:10:120)]) ; xlabel('frequency(hz)') ; ylabel('power(db)') ; 
subplot(2,2,3) ; imagesc(freqs,freqs,rstats,[-5,5]) ; title('t-tests for retinotopic tuning') ; colormap jet ; ylabel('frequency(hz)') ; xlabel('frequency(hz)') ; 
subplot(2,2,4) ; imagesc(freqs,freqs,ostats,[-5,5]) ; title('t-tests for orientation tuning') ; colormap jet ; ylabel('frequency(hz)') ; xlabel('frequency(hz)') ; 


figure,subplot(2,4,1);
shadedErrorBar(1:2:120,squeeze(mean(odprimes,1)),squeeze(std(odprimes,0,1))/3); xlim([0,120]); 
set(gcf,'Color','w'); xlabel('frequency(hz)') ; ylabel('orientation tuning'); 


rotvals = 1:400 ; 
shadedErrorBar([],squeeze(mean(mean(orvals(:,1,14:15,:),1),3)),squeeze(std(mean(orvals(:,1,14:15,:),3),0,1))./sqrt(5),{'m'}) ; hold on ; 
shadedErrorBar([],squeeze(mean(mean(orvals(:,1,4:8,:),1),3)),squeeze(std(mean(orvals(:,1,4:8,:),3),0,1))./sqrt(5),{'r'}) ; hline(0,'k') ; 
shadedErrorBar([],squeeze(mean(mean(orvals(:,1,26:30,:),1),3)),squeeze(std(mean(orvals(:,1,26:30,:),3),0,1))./sqrt(5),{'b'}) ; 
shadedErrorBar([],squeeze(mean(mean(orvals(:,1,10:12,:),1),3)),squeeze(std(mean(orvals(:,1,10:12,:),3),0,1))./sqrt(5),{'g'}) ; 
xlabel('orientation(deg)') ; ylabel('power(db)') ; set(gca,'XTick',1:45:size(orvals,4),'XTickLabel',[0,rotvals(45:45:end)]) ; 


clear allts 
for i=1:60 ; disp(i) ; 
    for j=1:26
        for k=1:26
            [h,p,ci,stats] = ttest2(squeeze(orvals(:,1,i,j)),squeeze(orvals(:,1,i,k))) ; 
            allts(i,j,k) = stats.tstat ; 
        end
    end
end




subplot(3,2,3) ; 
shadedErrorBar([],squeeze(mean(mean(retvals(:,1,8:10,:),1),3)),squeeze(std(mean(retvals(:,1,8:10,:),3),0,1))./sqrt(6),{'r'}) ; hold on ; 
shadedErrorBar([],squeeze(mean(mean(retvals(:,1,5:7,:),1),3)),squeeze(std(mean(retvals(:,1,5:7,:),3),0,1))./sqrt(6),{'b'}) ; hline(0,'k');
shadedErrorBar([],squeeze(mean(mean(retvals(:,1,28:30,:),1),3)),squeeze(std(mean(retvals(:,1,28:30,:),3),0,1))./sqrt(6),{'m'}) ; 
shadedErrorBar([],squeeze(mean(mean(retvals(:,1,14:16,:),1),3)),squeeze(std(mean(retvals(:,1,14:16,:),3),0,1))./sqrt(6),{'g'}); ylim([-5,5]) ; xlabel('angle(deg)') ; ylabel('power') ; 
subplot(3,2,4) ; 
shadedErrorBar([],squeeze(mean(mean(orvals(:,1,8:10,:),1),3)),squeeze(std(mean(orvals(:,1,8:10,:),3),0,1))./sqrt(6),{'r'}) ; hold on ; 
shadedErrorBar([],squeeze(mean(mean(orvals(:,1,5:7,:),1),3)),squeeze(std(mean(orvals(:,1,5:7,:),3),0,1))./sqrt(6),{'b'}) ; hline(0,'k');
shadedErrorBar([],squeeze(mean(mean(orvals(:,1,28:30,:),1),3)),squeeze(std(mean(orvals(:,1,28:30,:),3),0,1))./sqrt(6),{'m'}) ; 
shadedErrorBar([],squeeze(mean(mean(orvals(:,1,14:16,:),1),3)),squeeze(std(mean(orvals(:,1,14:16,:),3),0,1))./sqrt(6),{'g'}); ylim([-6,6]) ; xlabel('angle(deg)') ; ylabel('power') ; 
subplot(3,2,6) ; 
shadedErrorBar([],mean(odprimes),std(odprimes,0,1)./sqrt(5),{'k'}) ; ylim([0.4,1.4]) ; title('d-prime (orientation tuning)') ; set(gca,'XTick',1:10:60,'XTickLabel',[0,20:20:120]) ; xlabel('frequency(hz)') ; ylabel('tuning') ; hold on ; 
vline(5:7,'b') ; vline(8:10,'r') ; vline(28:30,'m') ; vline(14:16,'g') ;  shadedErrorBar([],mean(odprimes),std(odprimes,0,1)./sqrt(5),{'k'}) ; ylim([0.4,1.4]) ; title('d-prime (orientation tuning)') ; set(gca,'XTick',1:10:60,'XTickLabel',[0,20:20:120]) ; xlabel('frequency(hz)') ; ylabel('tuning') ; 
subplot(3,2,5) ; 
shadedErrorBar([],mean(rdprimes),std(rdprimes,0,1)./sqrt(5),{'k'}) ; ylim([0.4,1.4]) ; title('d-prime (retinotopic tuning)') ; set(gca,'XTick',1:10:60,'XTickLabel',[0,20:20:120]) ; xlabel('frequency(hz)') ; ylabel('tuning') ; hold on ; 
vline(5:7,'b') ; vline(8:10,'r') ; vline(28:30,'m') ; vline(14:16,'g') ;  shadedErrorBar([],mean(rdprimes),std(rdprimes,0,1)./sqrt(5),{'k'}) ; ylim([0.4,1.4]) ; title('d-prime (retinotopic tuning)') ; set(gca,'XTick',1:10:60,'XTickLabel',[0,20:20:120]) ; xlabel('frequency(hz)') ; ylabel('tuning') ; 
subplot(3,2,1) ; 
imagesc(1:400,freqs,squeeze(mean(retvals(1:end,:,:,:),1)),[-4,4]) ; colormap jet ; axis xy ; xlabel('angle(deg)') ; title('retinotopic mapping') ; ylabel('frequency(hz)') ; 
subplot(3,2,2) ; 
imagesc(1:400,freqs,squeeze(mean(orvals(1:end,:,:,:),1)),[-4,4]) ; colormap jet ; axis xy ; xlabel('angle(deg)') ; title('orientation tuning') ; ylabel('frequency(hz)') ; 

subplot(2,4,1)
imagesc(50:350,1:2:120,squeeze(mean(allors(:,:,50:350),1))) ; axis xy ; colormap jet ; colorbar; xlabel('orientation(deg'); ylabel('frequency(hz)'); title('grand average'); 
set(gcf,'Color','w'); 
xlabel('orientation (deg)'); ylabel('frequency(hz)'); 
subplot(2,4,2); shadedErrorBar(50:350,squeeze(mean(mean(allors(:,25:30,50:350),1),2)),squeeze(std(mean(allors(:,25:30,50:350),2),0,1))./3,{'Color','b'}); title('high gamma');
xlim([49,351])
subplot(2,4,3); shadedErrorBar(50:350,squeeze(mean(mean(allors(:,13:16,50:350),1),2)),squeeze(std(mean(allors(:,13:16,50:350),2),0,1))./3,{'Color','g'}); title('low gamma');
xlim([49,351])
subplot(2,4,4); shadedErrorBar(50:350,squeeze(mean(mean(allors(:,6:8,50:350),1),2)),squeeze(std(mean(allors(:,6:8,50:350),2),0,1))./3,{'r'}); title('alpha'); 
xlim([49,351])
bands = {6:8,13:16,25:30}; ors = 50:350; 
resors = allors(:,:,ors); 
for bnd=1:length(bands)
    bbnd = squeeze(mean(mean(resors(:,bands{bnd},:),1),2)); 
    maxind = find(bbnd==max(bbnd)); 
    minind = find(bbnd==min(bbnd)); 
    [h,ps(bnd),ci,stats] = ttest2(squeeze(mean(mean(resors(:,bands{bnd},maxind),2),3)),squeeze(mean(mean(resors(:,bands{bnd},minind),2),3)));
    ts(bnd) = stats.tstat; 
    
    for f=1:60
            [h,bandps(f,bnd),ci,stats] = ttest(squeeze(mean(mean(resors(:,f,maxind),2),3)),squeeze(mean(mean(resors(:,f,minind),2),3)));
            bandts(f,bnd) = stats.tstat;
    end
end
subplot(2,4,6); 
barwitherr(squeeze(std(mean(resors(:,bands{1},[maxind,minind]),2),0,1))/3,squeeze(mean(mean(resors(:,bands{1},[maxind,minind]),2),1)));
set(gca,'XTickLabel',{'highest','lowest'}); title(['p=',num2str(ps(1))]);
subplot(2,4,7);
barwitherr(squeeze(std(mean(resors(:,bands{2},[maxind,minind]),2),0,1))/3,squeeze(mean(mean(resors(:,bands{2},[maxind,minind]),2),1)));
set(gca,'XTickLabel',{'highest','lowest'}); title(['p=',num2str(ps(2))]);
subplot(2,4,8); 
barwitherr(squeeze(std(mean(resors(:,bands{3},[maxind,minind]),2),0,1))/3,squeeze(mean(mean(resors(:,bands{3},[maxind,minind]),2),1)));
set(gca,'XTickLabel',{'highest','lowest'}); title(['p=',num2str(ps(3))]);

degs = 90:45:315; 
figure,
subplot(1,2,1); 
plot(ors,squeeze(mean(resors(:,[13:16,25:30],:),2))'); hold on; 
shadedErrorBar(ors,squeeze(mean(mean(resors(:,[13:16,25:30],:),1),2)),squeeze(std(mean(resors(:,[13:16,25:30],:),2),0,1))./3) ; vline([90:45:315],'k'); hold on; 
for i=1:length(degs) ; text(degs(i),4.2,num2str(degs(i))); end ; xlabel('orientation(deg)'); ylabel('power(db)') ; 
set(gcf,'Color','w'); 
subplot(1,2,2); 
plot(ors,squeeze(mean(resors(:,[6:8],:),2))'); hold on; 
shadedErrorBar(ors,squeeze(mean(mean(resors(:,[6:8],:),1),2)),squeeze(std(mean(resors(:,[6:8],:),2),0,1))./3) ; vline([90:45:315],'k'); hold on; 
for i=1:length(degs) ; text(degs(i),4.2,num2str(degs(i))); end ; xlabel('orientation(deg)'); ylabel('power(db)') ; 
set(gcf,'Color','w'); 

figure,
for i=1:length(subs)
    subplot(2,5,i) ; 
    imagesc(5:350,1:2:120,squeeze(resors(i,:,:)),[-4,4]) ; 
    axis xy ; colormap jet; if i==1; xlabel('orientation(deg)'); ylabel('frequency(hz)');end
    title(['subject=',num2str(i)]);
end
set(gcf,'Color','w'); 

% do max-min analysis

clear bandps bandts
for f=1:60
    bbnd = squeeze(mean(mean(resors(:,f,:),1),2)); 
    maxind = find(bbnd==max(bbnd)); 
    minind = find(bbnd==min(bbnd)); 
    [h,bandps(f),ci,stats] = ttest(squeeze(mean(mean(resors(:,f,maxind),2),3)),squeeze(mean(mean(resors(:,f,minind),2),3)));
    bandts(f) = stats.tstat;
end
fs=1:2:120;
figure,subplot(2,4,1);
plot(1:2:120,(bandts)) ; xlim([1,120]) ; title('t-values'); 
subplot(2,4,2); 
plot(1:2:120,(bandps)) ; hline(0.05,'b');xlim([1,120]) ; title('p-values, *p<0.05'); 
for i=1:length(fs) ; if bandps(i) < 0.05 ; text(fs(i)-2,0.4,'*'); end; end


obliques = [135:90:315] - 50; 
cardinals = [90:90:270] - 50; 

meanob = squeeze(mean(resors(:,:,obliques),3)); 
meancard = squeeze(mean(resors(:,:,cardinals),3)); 

meanboth = zeros(10,60,2) ; meanboth(:,:,1) = meanob ; meanboth(:,:,2) = meancard; 

subplot(2,4,1); 
[h,p,ci,stats] = ttest(squeeze(mean(meanboth(:,[bands{1},bands{1}],1),2)),squeeze(mean(meanboth(:,[bands{1},bands{1}],2),2))); 
bar(squeeze(mean(meanboth(:,[bands{1},bands{1}],:),2))) ; xlabel('subjects') ; ylabel('power') ; legend({'oblique','cardinal'});xlim([0,11]);
subplot(2,4,2);
barwitherr(squeeze(std(mean(meanboth(:,[bands{1},bands{1}],:),2),0,1))/3,squeeze(mean(mean(meanboth(:,[bands{1},bands{1}],:),1),2))); 
title(['p=',num2str(p)]); set(gca,'XTickLabel',{'oblique','cardinal'});  ylabel('power(db)'); xlabel('orientation'); 
set(gcf,'Color','w'); 

subplot(2,4,3); 
[h,p,ci,stats] = ttest(squeeze(mean(meanboth(:,[bands{2},bands{3}],1),2)),squeeze(mean(meanboth(:,[bands{2},bands{3}],2),2))); 
bar(squeeze(mean(meanboth(:,[bands{2},bands{3}],:),2))) ; xlabel('subjects') ; ylabel('power') ; legend({'oblique','cardinal'}); xlim([0,11]);
subplot(2,4,4);
barwitherr(squeeze(std(mean(meanboth(:,[bands{2},bands{3}],:),2),0,1))/3,squeeze(mean(mean(meanboth(:,[bands{2},bands{3}],:),1),2))); 
title(['p=',num2str(p)]); set(gca,'XTickLabel',{'oblique','cardinal'});  ylabel('power(db)'); xlabel('orientation');
set(gcf,'Color','w'); 

figure,
subplot(1,3,1)
shadedErrorBar([],squeeze(mean(mean(retvals(:,1,8:10,:),1),3)),squeeze(std(mean(retvals(:,1,8:10,:),3),0,1))./sqrt(6),{'r'}) ; hold on ; title('beta'); xlabel('wedge position (deg)'); 
subplot(1,3,2);
shadedErrorBar([],squeeze(mean(mean(retvals(:,1,5:7,:),1),3)),squeeze(std(mean(retvals(:,1,5:7,:),3),0,1))./sqrt(6),{'b'}) ; hline(0,'k'); title('alpha')
subplot(1,3,3); 
shadedErrorBar([],squeeze(mean(mean(retvals(:,1,28:30,:),1),3)),squeeze(std(mean(retvals(:,1,28:30,:),3),0,1))./sqrt(6),{'m'}) ; title('high gamma'); 

figure,
subplot(1,3,1) ; 
imagesc(1:400,freqs,squeeze(mean(orvals(1:end,:,:,:),1)),[-5,5]) ; colormap jet ; 
axis xy ; xlabel('angle(deg)') ; title('orientation tuning') ; ylabel('frequency(hz)') ; 
subplot(1,3,2) ; 
imagesc(1:400,freqs,abs(squeeze(mean(orvals(1:end,:,:,:),1))),[0,5]) ; colormap jet ; 
axis xy ; xlabel('angle(deg)') ; title('ABS orientation tuning') ; ylabel('frequency(hz)') ; 

low_orr = squeeze(mean(mean(abs(orvals(:,1,(freqs>8 & freqs<25),:)),1),3)); 
high_orr = squeeze(mean(mean(abs(orvals(:,1,(freqs>30 & freqs<65),:)),1),3)); 
subplot(1,3,3); 
plot(low_orr(50:end-50),high_orr(50:end-50),'k.') ; lsline; xlabel('low freq'); ylabel('high freq'); ylim([1.8,2.7]) ; xlim([1.9,2.4]); title('ABS alpha/beta vs gamma tuning'); 

subplot(1,4,1); 
shadedErrorBar(freqs,mean(odprimes),std(odprimes,0,1)./sqrt(5),{'k'}) ; hold on ; xlim([0,120]); xlabel('frequency(hz)') ;ylabel('d-prime'); title('orientation d-prime');
subplot(1,4,2); 
shadedErrorBar(freqs,mean(rdprimes),std(rdprimes,0,1)./sqrt(5),{'k'}) ; hold on ; xlim([0,120]); xlabel('frequency(hz)') ;ylabel('d-prime'); title('retinotopy d-prime'); 
subplot(1,4,3) ; 
shadedErrorBar(freqs,mean(rdprimes),std(rdprimes,0,1)./sqrt(5),{'k'}) ; hold on ; xlim([0,120]); xlabel('frequency(hz)') ;ylabel('d-prime'); hold on; 
shadedErrorBar(freqs,mean(odprimes),std(odprimes,0,1)./sqrt(5),{'g'}) ;
subplot(1,4,4);
plot([1,1],'k') ; hold on ; plot([1,2],'g'); legend({'retinotopy','orientation'})


imagesc(1:400,freqs,squeeze(mean(retvals(1:end,:,:,:),1)),[-3,3]) ; colormap jet ; 
axis xy ; xlabel('angle(deg)') ; title('retinotopic mapping') ; ylabel('frequency(hz)') ; 



figure,
subplot(1,4,1) ; 
imagesc(1:400,freqs,squeeze(mean(retvals(1:end,:,:,:),1)),[-3,3]) ; colormap jet ;
axis xy ; xlabel('angle(deg)') ; title('retinotopic mapping') ; ylabel('frequency(hz)') ; 
subplot(1,4,2) ; 
imagesc(1:400,freqs,abs(squeeze(mean(retvals(1:end,:,:,:),1))),[0,3]) ; colormap jet ; 
axis xy ; xlabel('angle(deg)') ; title('ABS retinotopic mapping') ; ylabel('frequency(hz)') ; 

low_ret = squeeze(mean(mean(abs(retvals(:,1,(freqs>15 & freqs<25),:)),1),3)); 
high_ret = squeeze(mean(mean(abs(retvals(:,1,(freqs>30 & freqs<65),:)),1),3)); 
subplot(1,4,3); 
plot(low_ret(50:end-50),high_ret(50:end-50),'k.') ; lsline; xlabel('beta'); ylabel('gamma');ylim([.4,1.1]) ; xlim([.9,2]); title('ABS beta vs gamma tuning'); 

low_ret = squeeze(mean(mean(abs(retvals(:,1,(freqs>8 & freqs<13),:)),1),3)); 
high_ret = squeeze(mean(mean(abs(retvals(:,1,(freqs>30 & freqs<65),:)),1),3)); 
subplot(1,4,4); 
plot(low_ret(50:end-50),high_ret(50:end-50),'k.') ; lsline; xlabel('alpha'); ylabel('gamma'); ylim([.4,1.1]) ; xlim([1.4,3.5]); title('ABS alpha vs gamma tuning'); 


%}


%{
cd E:\orientation_retinotopy\sub_russell ; 
all_ors = load('all_ors'); all_ors = all_ors.all_ors;
or_axis = load('or_axis') ; or_axis = or_axis.or_axis; 

fors = imresize(squeeze(mean(mean(all_ors,1),2))',[1,360]); 

m_orvals=  squeeze(mean(orvals(:,:,:,1:360),1)); 

[fcorrs,fps] = corr(fors',m_orvals'); 

subplot(1,2,1); 
plot(squeeze(mean(m_orvals([8:12],:),1)),fors,'o'); lsline; xlabel('EEG power') ;ylabel('BOLD increase'); title('EEG Beta vs BOLD tuning');
subplot(1,2,2); 
plot(squeeze(mean(m_orvals([20:60],:),1)),fors,'o'); lsline; xlabel('EEG power') ;ylabel('BOLD increase'); title('EEG Gamma vs BOLD tuning'); 

shadedErrorBar([],squeeze(mean(mean(orvals(:,1,14:16,:),1),3)),squeeze(std(mean(orvals(:,1,14:16,:),3),0,1))./sqrt(6),{'k'}); xlabel('angle(deg)') ; ylabel('power') ; 


figure,
subplot(1,2,1);
imagesc(sv,freqs,abs(m_orvals)) ; axis xy ; colormap jet ;
subplot(1,2,2)
imagesc(sv,freqs,m_orvals); axis xy ; colormap jet; 

%}

