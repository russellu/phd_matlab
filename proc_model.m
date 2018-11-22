cd c:/shared/eeg_paper_3/ ; ls 
a = load('Gamma -  EEG traces forr all contrasts and correlations.txt') ; 

clear epochs 
scount = 1 ; 
for i=1:40
    for j=1:40
        epochs(i,j,:) = downsample(a((scount-1)*2048+1:scount*2048,4),4) ; 
        scount = scount + 1 ; 
    end
end
clear ersp 
for i=1:40 ;
    for j=1:40
        [ersp(i,j,:,:),itc,powbase,times,freqs] = newtimef(squeeze(epochs(i,j,:)),size(epochs,3),[0,2],256,0,...
            'plotersp','off','plotitc','off','baseline',NaN,'freqs',[1,120],'nfreqs',60,'winsize',64) ; 
        
    end
end

mersp = squeeze(mean(ersp,4)) ; 
smersp = squeeze(mean(mersp,3)) ; 


clvls = unique(a(:,1)) ; rlvls = unique(a(:,2)) ; 
figure,

mcontrast = mat2gray(squeeze(mean(mersp,2))) ; 
mrand = mat2gray(squeeze(mean(mersp,1))) ; 


labs1 = {'unpt','33%contrast','5%contrast'} ; labs2 = {'unpt','10%random','60%random'} ; 
subplot(2,2,1) ; 
plot(mcontrast([40,14,3],:)','LineWidth',3) ; set(gca,'XTick',1:5:length(freqs),'XTickLabel',round(freqs(1:5:end))) ; legend(labs1) ;xlabel('frequency(hz)') ; ylabel('AU'); 
subplot(2,2,2) ; 
plot(mrand([40,36,17],:)','LineWidth',3) ; set(gca,'XTick',1:5:length(freqs),'XTickLabel',round(freqs(1:5:end))) ; legend(labs2) ; xlabel('frequency(hz)') ; ylabel('AU'); 
subplot(2,2,3) ; 
imagesc(squeeze(sum(mersp,3))) ; ylabel('contrast') ; xlabel('randomization') ;
set(gca,'XTick',1:5:length(rlvls),'XTickLabel',floor(rlvls(length(rlvls):-5:1)*10)/10,'YTick',1:5:length(clvls),'YTickLabel',floor(clvls(1:5:length(clvls))*10)/10) ;


figure,
subplot(2,1,1) ; plot(mat2gray(mean(smersp,2)),'k','LineWidth',3) ; xlabel('contrast level') ; 
set(gca,'XTick',1:5:length(rlvls),'XTickLabel',floor(clvls(1:5:length(clvls))*10)/10) ; ylabel('AU'); 
subplot(2,1,2) ;plot(mat2gray(mean(smersp,1)),'k','LineWidth',3) ; xlabel('randomization level') ;
set(gca,'XTick',1:5:length(clvls),'XTickLabel',floor(clvls(length(clvls):-5:1)*10)/10) ; ylabel('AU') ; 






% raster plots
rast = load('Sample Raster - 1s sham+2s stim.txt') ; 








