cd c:/shared/fpap ; ls 
subs = {'sub_alex','sub_charest','sub_esteban','sub_fabio','sub_gab','sub_gabriella','sub_genevieve','sub_gina','sub_jeremie','sub_julie','sub_katrine','sub_lisa'...
        ,'sub_marc','sub_marie','sub_mathieu','sub_maxime','sub_mingham','sub_patricia','sub_po','sub_russell','sub_sunachakan','sub_vincent'} ;
fmri = load('basestims.mat') ; fmri = fmri.basestims ; 
eeg = load('eegsaves.mat') ; eeg = eeg.eegsaves ; 
freqs = load('freqs') ;freqs = freqs.freqs ; times = load('times') ; times = times.times ; 
times_s = (times/1000) ; 
stims = [2,3,1,5,6] ; 
meanf = squeeze(mean(fmri,3)) ; 
dispfreqs = find(freqs>7) ; 
freqinds = 57-[7,17,26,37,47] ; freqlabs = [20,40,60,80,100] ; 
timeinds = [43,101,159] ; timelabs = [0,1,2] ; 

figh = figure ; 
for s=1:length(stims) ; 
   subplot(1,5,s) ; imagesc(flipud(squeeze(mean(eeg(:,stims(s),dispfreqs,:),1))),[-2,2]) ; 
   set(gca,'XTick',timeinds,'XTickLabel',(timelabs),'YTick',fliplr(freqinds),'YTickLabel',fliplr(freqlabs)) ; vline([timeinds(1),timeinds(3)],'k') ;    
   if s==1 ; 
   ylabel('frequency(hz)') ; xlabel('time(s)') ; 
   end
end
set(figh, 'Position', [100, 100, 1000, 200]);

% plot the FMRI component mask
cd c:/shared/regf1 ; comps = dir('comps_*') ; clear allcomps
for comp=1:length(comps)
   c = load_untouch_nii(comps(comp).name) ; 
   allcomps(comp,:,:,:) = c.img ; 
   f = load_untouch_nii(['reg_',subs{comp},'.nii.gz']) ; 
   allfs(comp,:,:,:) = f.img ; 
end
meanf = squeeze(mean(allfs,1)) ; 
meancomp = squeeze(mean(allcomps,1)) ; 
mask = mat2gray(meancomp) ; 
figure,
offset = 12 ; 
for indx=1:7 ; 
    h = subplot(1,7,indx) ; 
    plotoverlayIntensity2D(squeeze(meanf(:,:,indx+offset)),squeeze(mask(:,:,indx+offset)),squeeze(meancomp(:,:,indx+offset))) ; 
    p = get(h, 'pos') ; 
    if indx>1 
        p(1) = prevleft + prevwidth ; 
        set(h, 'pos', p);
    end
    prevleft = p(1) ; 
    prevwidth = p(3) ; 
end
suptitle('grand average FMRI component weights'); 

meeg = squeeze(mean(eeg(:,:,:,times<2000 & times>0),4)) ; 
mfmri = squeeze(mean(fmri,3)) ; 

stims1 = [1,3,2] ; stims2 = [1,5,6] ; 
stimcolors1 = {[0,0,0],[1,0,0],[.5,0,0]} ; stimcolors2 = {[0,0,0],[0,1,0],[0,.5,0]} ; 
gfreqs = find(freqs>35) ; labs1 = {'unperturbed','33%contrast','5%contrast'} ; labs2 = {'unperturbed','10%randomized','60%randomized'} ;
freqinds = [3,13,23,33] ; freqlabs = [40,60,80,100] ;  timeinds = 1:10 ; timelabs = -2:2:18 ; 
%%% figure EEG and BOLD 
subplot(2,2,1) ; 
for i=1:length(stims1)
    errorbar(squeeze(mean(meeg(:,stims1(i),gfreqs),1)),squeeze(std(meeg(:,stims1(i),gfreqs),0,1))./sqrt(22),'Color',stimcolors1{i},'LineWidth',2) ; hold on ;
    xlim([1,length(gfreqs)]) ; ylabel('power(db)') ; set(gca,'XTick',freqinds,'XTickLabel',freqlabs) ; xlabel('frequency(hz)') ; hline(0,'k') ; 
    legend(labs1) ; 
end
subplot(2,2,2) ; 
for i=1:length(stims1)
    errorbar(squeeze(mean(meeg(:,stims2(i),gfreqs),1)),squeeze(std(meeg(:,stims2(i),gfreqs),0,1))./sqrt(22),'Color',stimcolors2{i},'LineWidth',2) ; hold on ;
    xlim([1,length(gfreqs)]) ; ylabel('power(db)') ; set(gca,'XTick',freqinds,'XTickLabel',freqlabs) ; xlabel('frequency(hz)') ; hline(0,'k') ; 
    legend(labs2) ;
end
subplot(2,2,3) ; 
for i=1:length(stims1)
    errorbar(squeeze(mean(mfmri(:,stims1(i),:),1)),squeeze(std(mfmri(:,stims1(i),:),0,1))./sqrt(22),'Color',stimcolors1{i},'LineWidth',2) ; hold on ;  
    xlim([1,10]) ; hline(0,'k') ; vline([2,3],'r') ; text(2,1.4,'on') ; text(3,1.4,'off') ;  legend(labs1) ; set(gca,'XTick',timeinds,'XTickLabel',timelabs) ; xlabel('time(s)') 
    ylabel('component weight(unitless)')
end
subplot(2,2,4) ; 
for i=1:length(stims1)
    errorbar(squeeze(mean(mfmri(:,stims2(i),:),1)),squeeze(std(mfmri(:,stims2(i),:),0,1))./sqrt(22),'Color',stimcolors2{i},'LineWidth',2) ; hold on ;  
    xlim([1,10]) ; hline(0,'k') ; vline([2,3],'r') ; text(2,1.4,'on') ; text(3,1.4,'off') ;  legend(labs2) ; set(gca,'XTick',timeinds,'XTickLabel',timelabs) ; xlabel('time(s)') 
    ylabel('component weight(unitless)')
end

stims3 = [2,3,1,5,6] ; nbg = find(freqs>60 & freqs<=80) ; bbg = find(freqs>80 & freqs<110) ; 

subplot(1,3,1) ; barwitherr(squeeze(std(mean(meeg(:,stims3,nbg),3),0,1))./sqrt(22),squeeze(mean(mean(meeg(:,stims3,nbg),1),3))) ; 
ylabel('power(db)') ; xlabel('stimulus type') ; title('NBG (60-80Hz)') ;
subplot(1,3,2) ; barwitherr(squeeze(std(mean(meeg(:,stims3,bbg),3),0,1))./sqrt(22),squeeze(mean(mean(meeg(:,stims3,bbg),1),3))) ; 
ylabel('power(db)') ; title('BBG (80-110Hz)') ; 
subplot(1,3,3) ; barwitherr(squeeze(std(mean(mfmri(:,stims3,[5,6]),3),0,1))./sqrt(22),squeeze(mean(mean(mfmri(:,stims3,[5,6]),1),3))) ; ylim([0.75,2]) ;
ylabel('component weight(unitless)') ; title('BOLD (6-8s post stimulus)') ; 

%%%% stats :
anova1(squeeze(mean(meeg(:,[1,2],nbg),3)),[],'off') 
anova1(squeeze(mean(mfmri(:,[1,2],[5,6]),3)),[],'off')
anova1(squeeze(mean(meeg(:,[1,6],nbg),3)),[],'off') 
anova1(squeeze(mean(mfmri(:,[1,6],[5,6]),3)),[],'off')
anova1(squeeze(mean(meeg(:,[5,6],nbg),3)),[],'off') 
anova1(squeeze(mean(mfmri(:,[5,6],[5,6]),3)),[],'off')

figure,
subplot(2,2,1) ; 
eeg1 = squeeze(mean(mean(meeg(:,stims3,nbg),1),3)) ; bold1 = squeeze(mean(mean(mfmri(:,stims3,5:6),1),3)) ;
p = polyfit(eeg1,bold1,1) ; f = polyval(p,eeg1) ;
errorbarxy(eeg1,bold1,squeeze(std(mean(meeg(:,stims3,nbg),3),0,1))./sqrt(22),squeeze(std(mean(mfmri(:,stims3,5:6),3),0,1))./sqrt(22),{'o','k','k'}) ; hold on ; 
plot(eeg1,f,'r') ; xlim([0,1.2]) ; xlabel('EEG(db)') ; ylabel('BOLD (arb. units)') ;
[pi,p,ol] = Shepherd(eeg1',bold1',10000) ; title(['Shepherds pi = ',num2str(pi),' p=',num2str(p)]) ;

subplot(2,2,2) ; 
eeg1 = squeeze(mean(mean(meeg(:,stims3,bbg),1),3)) ; bold1 = squeeze(mean(mean(mfmri(:,stims3,5:6),1),3)) ;
p = polyfit(eeg1,bold1,1) ;f = polyval(p,eeg1) ;
errorbarxy(squeeze(mean(mean(meeg(:,stims3,bbg),1),3)),squeeze(mean(mean(mfmri(:,stims3,5:6),1),3)),...
    squeeze(std(mean(meeg(:,stims3,bbg),3),0,1))./sqrt(22),squeeze(std(mean(mfmri(:,stims3,5:6),3),0,1))./sqrt(22),{'o','k','k'}) ; hold on ; 
plot(eeg1,f,'r') ;xlabel('EEG(db)') ; ylabel('BOLD (arb. units)') ;
[pi,p,ol] = Shepherd(eeg1',bold1',10000) ;  title(['Shepherds pi = ',num2str(pi),' p=',num2str(p)]) ;

for i=1:22 ; 
    for j=1:60
        subcorrs(i,j) = corr2(meeg(i,stims3,j),squeeze(mean(mfmri(i,stims3,5:6),3))) ; 
    end
end
subplot(
errorbar(squeeze(mean(subcorrs(:,dispfreqs),1)),squeeze(std(subcorrs(:,dispfreqs),0,1)./sqrt(22)),'LineWidth',2) ; hline(0,'k') ; xlim([0,58]) ; 
set(gca,'XTick',[7,17,27,37,47],'XTickLabel',[20,40,60,80,100]) ; ylabel('mean r-value') ; xlabel('frequency(hz)') ; 
