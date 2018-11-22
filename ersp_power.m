clear all ; close all  ;
cd c:/shared/resmerged ; 
subs=dir('*') ; subs(1:2) = [] ; 

goodcs = {[20,7,26,15],[13,3,16,23],[12,7,3,11,9],[13,19,6,14],[13,11],[26,22,17,8,42],[14,27,23],[17,16,7,4],[20,7,12],[13,10,11,21],[7,16],[16,11,9,7],...
    [8,12,28,39,50,51],[8,16,19,33,40],[20,22,27],[13,7,12,27],[7,3,2,42],[22,13,12],[14,5,15,17,21],[4,6,7,27],[2,6,7,12,23],[19,15,8,3],[15,7,17,8],...
    [7,8,14],[4,15,32,37],[6,5,12,28],[7,1,17,3],[22,14,28],[10,9,16,18,54],[7,3,2]} ; 
badts = {[6],[],[],[2],[2,39],[1],[1,90],[55],[],[128],[],[],[90],[],[],[119],[104],[],[],[],[],[17,21,83],[110],[],[],[],[37],[],[],[]} ; 

allcs = {[5,7,8,9,14,15,18,19,20,26],[3,5,7,8,10,11,13,16,23,25,26,31,34,38],[3,4,5,7,8,9,11,12,13,15,16,18,20,21,27],[6,7,8,10,11,12,13,14,16,18,19,29],...
    [9,11,13,16,25,52,53],[3,4,5,6,8,17,22,23,26,35,42,43,44,54],[6,8,10,12,13,14,20,23,27,30],[2,4,5,7,13,16,17,28,37],[7,10,12,15,16,17,18,20,22,25,27,29,31,32,35,38],...
    [3,5,7,10,11,13,15,21,22,23,32],[3,4,7,8,6,10,11,13,15,16,20,22,26],[5,7,9,11,12,15,16,23,31,38,54],[6,8,12,14,17,20,26,28],[1,3,4,5,6,8,9,12,16,18,19,20,31,33],...
    [2,7,10,11,20,27,28],[7,8,10,11,12,13,15,17,22,27,29,32],[1,2,3,5,7,8,12,16,18,21],[1,12,13,16,17,19,21,22,30,31],[3,5,6,7,8,9,10,11,12,13,14,15,17,20,21,22,34],...
    [3,4,5,6,7,8,13,17,21,26,27,28,40,41,42,46],[2,4,6,7,9,12,23,27,28,32,33,40,41],[3,5,10,11,15,17,19,20],[7,8,9,10,11,15,16,17,18,19,22,23,24,27,33],...
    [1,7,8,9,14,16],[4,5,6,7,10,11,12,13,15,17,20,21,32,33,37],[2,3,5,6,9,12,13,14,17,18,22,28,33,35],[1,3,4,6,7,11,14,16,17,18,30,37],...
    [3,4,5,6,8,9,14,15,16,20,22,28],[1,6,7,8,9,10,11,16,18,24],[1,2,3,6,7,8,11,13,14,16,17,21,23]};

for s=1:length(subs)
    cd(['c:/shared/resmerged/',subs(s).name]) ; ls 
    EEG = pop_loadset('merged.set') ; fulldat = EEG.data ; 
    saveica = load('saveica.mat') ; saveica = saveica.saveica ; 
    weights = saveica{1} ; sphere = saveica{2} ;
    winv = pinv(weights*sphere) ; 
    %unpt
    acts = weights*sphere*EEG.data ; 
    EEG.data = acts ; trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'} ; 
    for ss=1%:6
    epi = pop_epoch(EEG,{trigs{ss}},[-.65,3]) ; 
    clear ersp ; 
    for i=1:length(goodcs{s})
        for j=1:135
            [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(epi.data(goodcs{s}(i),:,j),epi.pnts,[epi.xmin,epi.xmax],epi.srate,0,...
               'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off','timesout',100) ; 
         %    [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(epi.data((i),:,j),epi.pnts,[epi.xmin,epi.xmax],epi.srate,0,...
         %       'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off','timesout',100) ; 
        end
    end   
    bersp = ersp - repmat(mean(ersp(:,:,:,times<0),4),[1,1,1,100]) ; 
    allbersp(s,ss,:,:,:) = squeeze(mean(ersp,1)) ; 
    bads = zeros(1,135) ; bads(badts{s}) = 1 ; goods = find(bads==0) ; 
    mersp = squeeze(mean(mean(bersp(:,goods,:,:),1),2)) ; 
    allmersp(s,ss,:,:) = mersp ; figure,imagesc(mersp) ; 
    end

    
  
    %{
    figure,for i=1:64 ; subplot(10,13,i) ; topoplot(winv(:,i),EEG.chanlocs,'electrodes','off') ; title(i) ; end ; 
    [a,f] = spectopo(EEG.data(:,5000:150000),0,256,'plot','off') ; 
    maxa = max(a,[],2) ;mina = min(a,[],2) ; 
    for i=1:64 ; subplot(10,13,i+65) ; plot(a(i,:)) ; title(i) ; xlim([1,length(f)]) ; ylim([mina(i),maxa(i)]) ; set(gca,'XTickLabel',[],'YTickLabel',[]) ; end ; suptitle(subs(s).name) ; 
    %}
    %figure,for i=1:135 ; subplot(10,14,i) ; imagesc(squeeze(mersp(i,:,:)),[-10,10]) ; title(i) ; set(gca,'XTickLabel',[],'YTickLabel',[]) ; end ; suptitle(subs(s).name) ; 

    
    %figure,for i=1:64 ; subplot(10,13,i) ; imagesc(squeeze(ersp(i,:,:)),[-5,5]) ; set(gca,'XTick',[],'YTick',[]) ; title(i) ; end ; 
    %for i=1:64 ; subplot(10,13,i+65) ; topoplot(winv(:,i),EEG.chanlocs,'electrodes','off') ; title(i) ; end ; suptitle(subs(s).name) ; 
    %allersp(s,:,:) = mean(ersp(goodcs{s},:,:)) ; 
    %figure,imagesc(squeeze(mean(ersp(goodcs{s},:,:),1)),[-2,2]) ; 

   
    comps = goodcs{s} ; zs = zeros(1,64) ; zs(comps) = 1 ; bads = find(zs==0) ; 
    w = weights*sphere ; 
    acts = weights*sphere*fulldat ; 
    acts(bads,:) = 0 ; 
    invacts = pinv(w)*acts ; 
    specacts = invacts - eegfiltfft(invacts,EEG.srate,59,61) ; 
    gammaspec = eegfiltfft(specacts,EEG.srate,0,128) ; 
    EEG.data = gammaspec ;     elecs = [56,24,57,25,58,26,59,61,62,63,60,64,29,30,31] ; 
    for t=1:length(trigs)
        gammaeps = pop_epoch(EEG,{trigs{t}},[-.8,2.8]) ;
        meanelecs(s,t,:,:) = squeeze(mean(gammaeps.data(elecs,:,1:135),1)) ; 
    end
    %{
    EEG.data = invacts ; 
    alleps = pop_epoch(EEG,{'S 11'},[-.8,2.8]) ; 
    gpow = abs(gammaeps.data) ; 
    gpow = gpow(:,:,goods) - repmat(mean(gpow(:,gammaeps.times<0,goods),2),[1,size(gpow,2),1]) ; 
    mgpow = squeeze(mean(mean(gpow(:,gammaeps.times>0 & gammaeps.times<2000,:),2),3)) ; 
    figure,
    topoplot(mgpow,EEG.chanlocs) ; 
    elecs = [56,24,57,25,58,26,59,61,62,63,60,64,29,30,31] ; 
    subeps(s,:,:) = squeeze(mean(alleps.data(elecs,:,:),1)) ; 
    alltopos(s,:) = mgpow ; 
    %}
end

mbersp = allbersp - repmat(mean(allbersp(:,:,:,:,times<0),5),[1,1,1,1,100]) ; 

shadedErrorBar([],squeeze(mean(mean(mbersp(2,1,:,:,times>0 & times<2),3),5)),squeeze(std(mean(mbersp(2,1,:,:,times>0 & times<2),5),0,3))/sqrt(135),{'r'}) ; hold on ; 
shadedErrorBar([],squeeze(mean(mean(mbersp(1,1,:,:,times>0 & times<2),3),5)),squeeze(std(mean(mbersp(1,1,:,:,times>0 & times<2),5),0,3))/sqrt(135),{'b'}) ; 
shadedErrorBar([],squeeze(mean(mean(mbersp(7,1,:,:,times>0 & times<2),3),5)),squeeze(std(mean(mbersp(7,1,:,:,times>0 & times<2),5),0,3))/sqrt(135),{'g'}) ; 
set(gca,'XTick',1:5:length(freqs),'XTickLabel',round(freqs(1:5:end))) ; xlabel('frequency(hz)') ; ylabel('power(db)') ; 
figure,plot(1,'r') ; hold on ; plot(2,'b') ; plot(3,'g') ; legend({'subject 1','subject 2','subject 3'}) ; 


subplot(2,1,1) ; plot(epi.data(28,:,1)) ; subplot(2,2,3) ; imagesc(squeeze(bersp(28,1,:,:)),[-15,15]) ; axis xy ;

stinds = find(epi.times<0) ; endinds = find(epi.times<2000) ; 
subplot(2,1,1) ; plot(epi.times,epi.data(28,:,1)) ; vline([0,2000],'k') %set(gca,'XTick',1:50:length(epi.times),'XTickLabel',round(epi.times(1:50:end))) ; 
xlabel('time(ms)')  ;ylabel('signal(a.u)') ; xlim([epi.times(1),epi.times(end)])
subplot(2,1,2) ; imagesc(times,freqs,squeeze(bersp(28,1,:,:)),[-15,15]) ; axis xy ; xlabel('time(s)') ; ylabel('frequency(hz)') ; vline([0,2],'k')

baselecs = meanelecs - repmat(mean(meanelecs(:,:,1:80,:),3),[1,1,922,1]) ; 
melecs = squeeze(mean(baselecs,4)) ; 

shadedErrorBar([],squeeze(mean(melecs(:,1,:),1)),squeeze(std(melecs(:,1,:),0,1))/sqrt(30),{'r'}) ; hold on ; 
shadedErrorBar([],squeeze(mean(melecs(:,3,:),1)),squeeze(std(melecs(:,3,:),0,1))/sqrt(30),{'g'}) ; 
shadedErrorBar([],squeeze(mean(melecs(:,2,:),1)),squeeze(std(melecs(:,2,:),0,1))/sqrt(30),{'b'}) ; ylabel('microvolt') ; 
set(gca,'XTick',1:60:length(gammaeps.times),'XTickLabel',round(gammaeps.times(1:60:end))) ; xlabel('time(ms)') ; 

mmersp = squeeze(mean(mean(allmersp(:,:,:,times>0 & times<2),2),4)) ; 
mmelecs = squeeze(mean(melecs,2)) ; 
corrs = corr(mmersp,mmelecs) ; 
cmersp = corr(mmersp)
imagesc(freqs,freqs,cmersp,[-1,1]) ; xlabel('frequency(hz)') ; ylabel('frequency(hz)') ; 
correrp = corr(mmersp,mean(mmelecs(:,200:600),2)) ; 

submat = zeros(5,6) ; 
for i=1:30 ; 
    subplottight(6,5,i) ; imagesc(times(1:90),freqs,squeeze(mean(allmersp(i,[1,3,5],:,1:90),2)),[-3,3]) ; 
    axis xy ; set(gca,'XTick',[],'YTick',[]) ; text(2,110,['#',num2str(i)]) ; 
    [x,y] = ind2sub(size(submat),i) ; 
    submat(x,y) = squeeze(mean(mean(mean(allmersp(i,[1,3,5],freqs>50 & freqs<75,1:90),2),3),4)) ;
end


plot(squeeze(std(mean(mean(allmersp(:,1,:,times>0 & times<2),2),4),0,1)))
set(gca,'XTick',1:5:60,'XTickLabel',round(freqs(1:5:60))) ; xlabel('frequency(hz)') ; ylabel('std(db)') ; 


%subeps = squeeze(subeps) ; 
%{
baseps = subeps - repmat(mean(subeps(:,1:200,:),2),[1,922,1]) ; 
mtopos = squeeze(mean(alltopos(:,elecs),2)) ; 
for i=1:31
    bads = zeros(1,135) ; bads(badts{i}) = 1 ; goods = find(bads==0) ;
    meps(i,:) = mean(baseps(i,:,goods),3) ; 
end

mineps = min(meps(:,235:260),[],2) ; maxeps = max(meps(:,245:280),[],2) ; 
deps = maxeps - mineps ; clear corrs
for i=1:60
    for j=1:200
        corrs(i,j) = corr2(deps,allmersp(:,i,j)) ; 
    end
end
%}
% check single trial adn other stimulus types (across stimulus types, and
% trials) 
%{
gender = [0,1,1,0,0,0,0,1,1,1,0,0,1,1,1,1,0,1,0,0,1,0,1,0,0,1,1,1,1,0,1] ; 

[sv,si] = sort(mean(mean(allmersp(:,25:50,times>0 & times<2),2),3),'descend') ; 
newmersp = allmersp(si,:,:) ; gender = gender(si) ; 

allmale = squeeze(mean(newmersp(gender==0,:,:),1)) ; 
allfem = squeeze(mean(newmersp(gender==1,:,:),1)) ; 

subplot(1,2,1) ; imagesc(times,freqs,allmale,[-3,3]) ; axis xy ; 
subplot(1,2,2) ; imagesc(times,freqs,allfem,[-3,3]) ; axis xy ; set(gcf,'color','w') ; ylabel('frequency(hz)') ; xlabel('time(s)') ; 

mans = newmersp(find(gender==0),:,:) ; 
fems = newmersp(find(gender==1),:,:) ; 
for i=1:14 ; subplot(2,17,i) ; imagesc(times,freqs,squeeze(mans(i,:,:)),[-4,4]) ; axis xy ; end
for i=1:17 ; subplot(2,17,i+17) ; imagesc(times,freqs,squeeze(fems(i,:,:)),[-4,4]) ; axis xy ; end

for i=1:1:60
    for j=1:200
        [h,p,ci,stats] = ttest2(squeeze(newmersp(gender==0,i,j)),squeeze(newmersp(gender==1,i,j))) ; 
        allts(i,j) = stats.tstat ; 
    end
end

figure,subplot(2,4,1) ; 
imagesc(times,freqs,squeeze(std(mans,0,1)),[-2,2]) ; axis xy ; 
subplot(2,4,2) ; 
imagesc(times,freqs,squeeze(std(fems,0,1)),[-2,2]) ; axis xy ; 
subplot(2,4,3) ; 
imagesc(times,freqs,abs(squeeze(mean(fems,1))-squeeze(mean(mans,1))),[-2,2]) ; axis xy ; 
subplot(2,4,4) ; 
imagesc(times,freqs,allts,[-3,3]) ; axis xy ; 
%}
%{

subplot(3,3,1) ; imagesc(times,freqs,squeeze(mean(allmersp(:,1,:,:),1)),[-2,2]) ; axis xy ; xlabel('time(s)') ; ylabel('freq(hz)') ; title('100% contrast') ; 
subplot(3,3,2) ; imagesc(times,freqs,squeeze(mean(allmersp(:,2,:,:),1)),[-2,2]) ; axis xy ; title('33% contrast') ; 
subplot(3,3,3) ; imagesc(times,freqs,squeeze(mean(allmersp(:,3,:,:),1)),[-2,2]) ; axis xy ; title('5% contrast') ; 
subplot(2,1,2) ; shadedErrorBar([],mean(mean(allmersp(:,1,:,times>0 & times<2),1),4),std(mean(allmersp(:,1,:,times>0 & times<2),4),0,1)/sqrt(31),{'r'}) ; hold on ; 
shadedErrorBar([],mean(mean(allmersp(:,2,:,times>0 & times<2),1),4),std(mean(allmersp(:,2,:,times>0 & times<2),4),0,1)/sqrt(31),{'g'}) ; 
shadedErrorBar([],mean(mean(allmersp(:,3,:,times>0 & times<2),1),4),std(mean(allmersp(:,3,:,times>0 & times<2),4),0,1)/sqrt(31),{'b'}) ; 
set(gca,'XTick',1:5:60,'XTickLabel',round(freqs(1:5:60))) ; xlabel('freq(hz)') ;ylabel('power(db') ; 
figure,plot([1:3],'r') ; hold on ; plot([2:4],'g') ; plot(3:5,'b') ; legend({'100%','33%','5%'}) ; 
figure,imagesc([-2,2]) ; 


subplot(3,3,1) ; imagesc(times,freqs,squeeze(mean(allmersp(2,1,:,:),1)),[-4,4]) ; axis xy ; xlabel('time(s)') ; ylabel('freq(hz)') ; title('100% contrast') ; 
subplot(3,3,2) ; imagesc(times,freqs,squeeze(mean(allmersp(2,2,:,:),1)),[-4,4]) ; axis xy ; title('33% contrast') ; 
subplot(3,3,3) ; imagesc(times,freqs,squeeze(mean(allmersp(2,3,:,:),1)),[-4,4]) ; axis xy ; title('5% contrast') ; 

subplot(3,3,4) ; imagesc(times,freqs,squeeze(mean(allmersp(13,1,:,:),1)),[-4,4]) ; axis xy ; 
subplot(3,3,5) ; imagesc(times,freqs,squeeze(mean(allmersp(13,2,:,:),1)),[-4,4]) ; axis xy ;
subplot(3,3,6) ; imagesc(times,freqs,squeeze(mean(allmersp(13,3,:,:),1)),[-4,4]) ; axis xy ; 

subplot(3,3,7) ; imagesc(times,freqs,squeeze(mean(allmersp(1,1,:,:),1)),[-4,4]) ; axis xy ; 
subplot(3,3,8) ; imagesc(times,freqs,squeeze(mean(allmersp(1,2,:,:),1)),[-4,4]) ; axis xy ;
subplot(3,3,9) ; imagesc(times,freqs,squeeze(mean(allmersp(1,3,:,:),1)),[-4,4]) ; axis xy ; 
%}

%{
mtersp = squeeze(mean(mean(allmersp(:,[1,2,3],:,times>1 & times<2),2),4)) ; 
[c,p] = corr(mtersp) ; 
stim1 = (squeeze(mean(allmersp(:,1,:,:),2))) ; 
for i=1:200 ; corrs(i,:,:)  =corr(squeeze(stim1(:,:,i))) ; end
figure,for i=1:200 ; subplottight(10,20,i) ; imagesc(squeeze(corrs(i,:,:)),[-.75,.75]) ; end
subplot(1,2,1) ; imagesc(squeeze(mean(corrs(times>0 & times<0.15,:,:),1)),[-.8,.8]) ; 
subplot(1,2,2) ; imagesc(squeeze(mean(corrs(times>1.5 & times<2,:,:),1)),[-.8,.8]) ; 
figure,
plot(squeeze(corrs(:,3,28))) ; hold on ; hline(0,'k') ; plot(corrs(:,54,46),'r') ; 
%}

tt = times>2.5 & times<3 ; 
subplot(2,2,1) ; 
shadedErrorBar([],mean(mean(allmersp(:,1,:,tt),1),4),std(mean(allmersp(:,1,:,tt),4),0,1)/sqrt(31),{'r'}) ; hold on ; 
shadedErrorBar([],mean(mean(allmersp(:,3,:,tt),1),4),std(mean(allmersp(:,3,:,tt),4),0,1)/sqrt(31),{'g'}) ; 
shadedErrorBar([],mean(mean(allmersp(:,2,:,tt),1),4),std(mean(allmersp(:,2,:,tt),4),0,1)/sqrt(31),{'b'}) ;  hline(0,'k') ; 
subplot(2,2,2) ; 
shadedErrorBar([],mean(mean(allmersp(:,1,:,tt),1),4),std(mean(allmersp(:,1,:,tt),4),0,1)/sqrt(31),{'r'}) ; hold on ; 
shadedErrorBar([],mean(mean(allmersp(:,5,:,tt),1),4),std(mean(allmersp(:,5,:,tt),4),0,1)/sqrt(31),{'m'}) ; 
shadedErrorBar([],mean(mean(allmersp(:,6,:,tt),1),4),std(mean(allmersp(:,6,:,tt),4),0,1)/sqrt(31),{'c'}) ;  hline(0,'k') ; 
subplot(2,2,3) ; 
shadedErrorBar([],mean(mean(allmersp(:,1,:,tt),1),4),std(mean(allmersp(:,1,:,tt),4),0,1)/sqrt(31),{'r'}) ; hold on ; 
shadedErrorBar([],mean(mean(allmersp(:,4,:,tt),1),4),std(mean(allmersp(:,4,:,tt),4),0,1)/sqrt(31),{'k'}) ; hline(0,'k') ; 

corrbersp = allbersp - repmat(mean(allbersp(:,:,:,:,times<0),5),[1,1,1,1,100]) ; 

for i=1:30 ; 
    figure('units','normalized','outerposition',[0 0 1 1])
    for j=1:135 ; subplot(10,14,j) ; 
        imagesc(squeeze(corrbersp(i,3,j,:,:)),[-10,10]) ; title(j); set(gca,'XTick',[],'YTick',[]) ; 
    end ; suptitle(subs(i).name) ;  
end ; 

bad1s = {[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]} ;
bad1s{28} = [80] ; bad1s{27} = [126,54,37,58] ; bad1s{26} = [123] ; bad1s{23} = [110,40] ; bad1s{22} = [17,83] ; bad1s{21} = [118] ; 
bad1s{20} = [34,116,127] ; bad1s{18} = [50,54] ; bad1s{17} = [104,47] ; bad1s{16} = [119,112,58] ; bad1s{13} = [34] ; bad1s{12} = [98] ; 
bad1s{10} = [128,103] ; bad1s{8} = [55] ; bad1s{5} = [39] ; bad1s{4} = [2] ; bad1s{3} = [54,73,118] ; bad1s{1} = [6] ; 
bad2s = {[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]} ;
bad2s{30} = [62] ; bad2s{29} = [44] ; bad2s{28} = [83] ; bad2s{26} = [1,16,105] ; bad2s{23} = 50 ; bad2s{20} = [65,93] ; bad2s{19} = [18,22,42,132] ; 
bad2s{17} = [77,68] ; bad2s{13} = [67,121] ; bad2s{12} = [106,64] ; bad2s{9} = [106,22] ; bad2s{4} = [2,57] ; bad2s{1} = [29,65] ;  

bad3s = {[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]} ;
bad4s = {[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]} ;
bad5s = {[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]} ;
bad6s = {[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]} ;


