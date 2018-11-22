clear all ; close all  ;
cd c:/shared/papsaves ; times = load('times') ; times = times.times ; freqs = load('freqs') ; freqs = freqs.freqs ; 

% fmri 
allstimepochs = load('allstimepochs') ; allstimepochs = allstimepochs.allstimepochs ; 
for i=1:22 
    for j=1:6 ;
        for k=1:45
            allstimepochs(i,j,k,:) = allstimepochs(i,j,k,:) - squeeze(allstimepochs(i,j,k,2)) ; 
        end
    end
end

mstims = squeeze(mean(mean(allstimepochs(:,stims,:,5:6),3),4)) ; 

%%% EEG
allersp = load('allersp.mat') ; allersp = allersp.allersp ; 
times = load('times') ; times = times.times ; freqs = load('freqs') ; freqs = freqs.freqs ; 
stims = [2,3,1,5,6] ;
%%% baseline correct
%{
for i=1:22 ; 
    for j=1:6 ; 
        for k=1:6 ;
            for el=1:135;
                basijkl = repmat(squeeze(mean(allersp(i,j,k,el,:,times<0),6)),[1,200]) ; 
                allersp(i,j,k,el,:,:) = (squeeze(allersp(i,j,k,el,:,:))-basijkl) ; 
            end
        end
    end
end
%}

allersp = allersp - repmat(mean(allersp(:,:,:,:,:,times<0),6),[1,1,1,1,1,200]) ; 

subj = squeeze(allersp(19,:,:,:,:,:)) ; 


% remove bad trials
mersp = squeeze(mean(mean(mean(allersp(:,:,1:4,:,freqs>5 & freqs<110,times>0 & times<2000),3),5),6)) ;

mersp2 = squeeze(mean(mean(allersp(:,:,1:4,:,:,:),3),4)) ; 
stmersp = mtersp(:,stims,:,:) ; 

for i=1:22 ; 
    for j=1:6
       zsij = abs(zscore(squeeze(mersp(i,j,:)))) ; 
       bads(i,j,:) = zsij>=2 ; 
    end
end

% get the mean ersp across trials
mcersp = squeeze(mean(allersp(:,:,1:4,:,:,:),3)) ; clear mtersp
for i=1:22
    for j=1:6
        mtersp(i,j,:,:) = squeeze(mean(mcersp(i,j,bads(i,j,:)==0,:,:),3)) ; 
    end
end


%%% get significant differences across stimuls types for all EEG subjects
nbg = find(freqs>50 & freqs<80) ; t = find(times>0 & times<2) ; 
tstmersp = squeeze(mean(mean(mean(allersp(:,:,1:2,:,nbg,t),3),5),6)) ; 
stimtypes = [2,3,4,5,6] ; 
for i=1:size(allersp,1)
    for j=1:length(stimtypes) ; 
        [h(i,j),p(i,j),ci,stats] = ttest2(squeeze(tstmersp(i,1,squeeze(bads(i,1,:))==0)),squeeze(tstmersp(i,stimtypes(j),squeeze(bads(i,stimtypes(j),:))==0))) ; 
        ts(i,j) = stats.tstat ; 
    end
end
 figure,
subplot(1,4,1) ;bar(ts(:,1)) ; set(gca,'XTick',1:size(h,1),'XTickLabel',h(:,1)) ; xlabel('subjects, 0/1 => p>0.05/p<0.05') ; ylabel('t-value (indep. samples)') ; xlim([0,23]) ; 
title(['unperturbed-5%contrast p-values: ',num2str(sum(h(:,1))),'/',num2str(size(h,1)),' subjects significant']) ; 
subplot(1,4,2) ;bar(ts(:,2)) ; set(gca,'XTick',1:size(h,1),'XTickLabel',h(:,2)) ; xlabel('subjects, 0/1 => p>0.05/p<0.05') ; ylabel('t-value (indep. samples)') ; xlim([0,23]) ;
title(['unperturbed-33%contrast p-values: ',num2str(sum(h(:,2))),'/',num2str(size(h,1)),' subjects significant']) ; 
subplot(1,4,3) ;bar(ts(:,4)) ; set(gca,'XTick',1:size(h,1),'XTickLabel',h(:,4)) ; xlabel('subjects, 0/1 => p>0.05/p<0.05') ; ylabel('t-value (indep. samples)') ; xlim([0,23]) ;
title(['unperturbed-10%random p-values: ',num2str(sum(h(:,4))),'/',num2str(size(h,1)),' subjects significant']) ; 
subplot(1,4,4) ;bar(ts(:,5)) ; set(gca,'XTick',1:size(h,1),'XTickLabel',h(:,5)) ; xlabel('subjects, 0/1 => p>0.05/p<0.05') ; ylabel('t-value (indep. samples)') ; xlim([0,23]) ;
title(['unperturbed-60%random p-values: ',num2str(sum(h(:,5))),'/',num2str(size(h,1)),' subjects significant']) ; 


%%% bold t-values
mstimepochs = squeeze(mean(allstimepochs(:,:,:,5:6),4)) ; 
stimtypes = [2,3,4,5,6] ; 
for i=1:size(allstimepochs,1)
    for j=1:length(stimtypes)
        [h(i,j),p(i,j),ci,stats] = ttest2(squeeze(mstimepochs(i,1,:)),squeeze(mstimepochs(i,stimtypes(j),:))) ; 
        ts(i,j) = stats.tstat ; 
    end
end
subplot(1,4,1) ;bar(ts(:,1)) ; set(gca,'XTick',1:size(h,1),'XTickLabel',h(:,1)) ; xlabel('subjects, 0/1 => p>0.05/p<0.05') ; ylabel('t-value (indep. samples)') ; xlim([0,23]) ; 
title(['unperturbed-5%contrast p-values: ',num2str(sum(h(:,1))),'/',num2str(size(h,1)),' subjects significant']) ; 
subplot(1,4,2) ;bar(ts(:,2)) ; set(gca,'XTick',1:size(h,1),'XTickLabel',h(:,2)) ; xlabel('subjects, 0/1 => p>0.05/p<0.05') ; ylabel('t-value (indep. samples)') ; xlim([0,23]) ;
title(['unperturbed-33%contrast p-values: ',num2str(sum(h(:,2))),'/',num2str(size(h,1)),' subjects significant']) ; 
subplot(1,4,3) ;bar(ts(:,4)) ; set(gca,'XTick',1:size(h,1),'XTickLabel',h(:,4)) ; xlabel('subjects, 0/1 => p>0.05/p<0.05') ; ylabel('t-value (indep. samples)') ; xlim([0,23]) ;
title(['unperturbed-10%random p-values: ',num2str(sum(h(:,4))),'/',num2str(size(h,1)),' subjects significant']) ; 
subplot(1,4,4) ;bar(ts(:,5)) ; set(gca,'XTick',1:size(h,1),'XTickLabel',h(:,5)) ; xlabel('subjects, 0/1 => p>0.05/p<0.05') ; ylabel('t-value (indep. samples)') ; xlim([0,23]) ;
title(['unperturbed-60%random p-values: ',num2str(sum(h(:,5))),'/',num2str(size(h,1)),' subjects significant']) ; 





% do stuff to compare the two
meanf = squeeze(mean(mean(allstimepochs(:,:,:,5:6),3),4)) ; 
stims = [2,3,1,5,6] ; 
for i=1:size(mtersp,3); 
    for j=1:size(mtersp,4) 
        corrersp(i,j) = corr2(squeeze(meanf(:,stims)),mtersp(:,stims,i,j)) ; 
        
    end
end


for i=1:size(mtersp,3); 
    for j=1:size(mtersp,4) 
        mcorrersp(i,j) = corr2(squeeze(mean(meanf(:,stims),1)),mean(mtersp(:,stims,i,j),1)) ; 
        
    end
end


for i=1:22 ; 
    for j=1:60
        for k=1:200
            cersp(i,j,k) = corr2(squeeze(mtersp(i,stims,j,k)),squeeze(meanf(i,stims))) ; 
            
        end
    end
end
gfreqs = find(freqs>8) ; trimfreqs = freqs(gfreqs) ; 
fbands = [20,40,60,80,100] ; freqinds = zeros(size(fbands)) ; 
for i=1:length(fbands) ; diffs = abs(freqs(gfreqs)-fbands(i)) ; freqinds(i) = find(diffs==min(diffs)) ; end
tbands = [0,1,2] ; tinds = [44,101,159] ; 
subplot(2,2,1) ; mcersp = squeeze(mean(cersp,1)) ;  
imagesc(mcersp(gfreqs,:)) ; 
set(gca,'YTick',freqinds,'YTickLabel',fbands,'XTick',tinds,'XTickLabel',tbands) ; axis xy
xlabel('time(s)') ; ylabel('frequency(hz)') ; c=colorbar ; title(c,'r') ;  vline([tinds(1),tinds(3)],'k') ;  

%%% get the single subjet BOLD from the correlation ERSP
mcersp = squeeze(mean(cersp,1)) ; 
for i=1:22 ;
    for j=1:6
        mbolds(i,j) = sum(sum(squeeze(mtersp(i,j,:,:)).*mcersp)) ; 
    end
end
% modeled BOLD
subplot(2,2,1) ; barwitherr(squeeze(std(mbolds(:,stims),0,1))./sqrt(22),squeeze(mean(mbolds(:,stims),1))) ; ylabel('modeled BOLD') ; 



gfreqs = find(freqs>8) ; trimfreqs = freqs(gfreqs) ; 
fbands = [20,40,60,80,100] ; freqinds = zeros(size(fbands)) ; 
for i=1:length(fbands) ; diffs = abs(freqs(gfreqs)-fbands(i)) ; freqinds(i) = find(diffs==min(diffs)) ; end
tbands = [0,1,2] ; tinds = [44,101,159] ; 

subplot(2,2,1) ; 
imagesc(mcorrersp(gfreqs,:)) ; 
set(gca,'YTick',freqinds,'YTickLabel',fbands,'XTick',tinds,'XTickLabel',tbands) ; axis xy
xlabel('time(s)') ; ylabel('frequency(hz)') ; c=colorbar ; title(c,'r') ;  vline([tinds(1),tinds(3)],'k') ;  

subs = {'sub_alex','sub_charest','sub_esteban','sub_fabio','sub_gab','sub_gabriella','sub_genevieve','sub_gina','sub_jeremie','sub_julie','sub_katrine','sub_lisa'...
        ,'sub_marc','sub_marie','sub_mathieu','sub_maxime','sub_mingham','sub_patricia','sub_po','sub_russell','sub_sunachakan','sub_vincent'} ;
%%% get the 5 stimulus types in MNI space
    for s=1:length(subs) ; 
       cd(['c:/shared/allfmris/',subs{s},'/trigs']) ;  
       t1vox = dir('mni_*') ; 
       disp(size(t1vox)) ; 
       disp(subs{s}); 
       for t1=1:length(t1vox) ; 
           curr = load_untouch_nii(t1vox(t1).name) ;
           allvox(s,t1,:,:,:) = curr.img ; 
       end      
    end

    
    
    corrvox = zeros(182,218,182)  ;
    meanbeta = squeeze(mean(mean(mean(mtersp(:,:,freqs>8 & freqs<25,times>0 & times<2),1),3),4)) ;
    mvox = squeeze(mean(allvox,1)) ; 
    for i=1:182 ; disp(i) ; 
        for j=1:218
            for k=1:182
                corrvox(i,j,k) = corr2(meanbeta(:,stims),squeeze(mvox(stims,i,j,k))') ;               
            end
        end
    end
    
    mni = load_untouch_nii('c:/shared/ATLASES/ss.nii.gz') ;  
    mnim = mni.img ;
    corrvox = corrvox.^8 ; 
    
    
   cvox = corrvox ; 
    scount = 1 ; 
    for slices = 60%:2:100 ; 
        subplot(1,1,scount)  ;
        scount = scount + 1 ; 
        plotoverlayIntensity2D(squeeze(mnim(:,:,slices)),squeeze(abs(cvox(:,:,slices))),squeeze(cvox(:,:,slices)),90) ; 
    end
    
    
    
    mmtersp = squeeze(mean(mtersp(:,:,:,times>0&times<2),4)) ; 
    icount = 1 ; 
    for i=1:2:60-5
        subplot(3,8,icount) ; icount=icount+1,plot(squeeze(mean(meanf(:,stims),1)),squeeze(mean(mean(mmtersp(:,stims,i:i+5),1),3)),'o') ; 
        title(num2str(corr2(squeeze(mean(meanf(:,stims),1)),squeeze(mean(mean(mmtersp(:,stims,i:i+5),1),3)))))
    end

    
    mmtersp = squeeze(mean(mtersp(:,:,:,times>0&times<2),4)) ; 
    icount = 1 ; 
    for i=1:5:60-5
        subplot(3,4,icount) ; icount=icount+1,plot(squeeze(meanf(:,stims)),squeeze(mean(mmtersp(:,stims,i:i+5),3)),'o') ; 
        title(num2str(corr2(squeeze(meanf(:,stims)),squeeze(mean(mmtersp(:,stims,i:i+5),3)))))
    end
    
    subplot(1,3,1) ; 
    hz = (freqs>8 & freqs<25) ; 
    x = squeeze(mean(mean(mmtersp(:,stims,hz),3),1)); 
    y = squeeze(mean(meanf(:,stims),1)); 
    P = polyfit(x,y,1);yfit = P(1)*x+P(2);   
    hz = find(freqs>80 & freqs<100 ) ;  plot(x,yfit,'r','LineWidth',2);hold on ; 
    errorbarxy(x,y,squeeze(std(mean(mmtersp(:,stims,hz),3),0,1))./sqrt(22),squeeze(std(meanf(:,stims),0,1)./sqrt(22)),{'o','k','k'}) ; 
    xlabel('alpha/beta (db)') ; ylabel('BOLD (arb.units)');
    subplot(1,3,2) ;
    hz = (freqs>55 & freqs<75) ; 
    x = squeeze(mean(mean(mmtersp(:,stims,hz),3),1)); 
    y = squeeze(mean(meanf(:,stims),1)); 
    P = polyfit(x,y,1);yfit = P(1)*x+P(2);   
    hz = find(freqs>80 & freqs<100 ) ;  plot(x,yfit,'r','LineWidth',2);hold on ; 
    errorbarxy(x,y,squeeze(std(mean(mmtersp(:,stims,hz),3),0,1))./sqrt(22),squeeze(std(meanf(:,stims),0,1)./sqrt(22)),{'o','k','k'}) ; 
    xlabel('NBG (db)') ; ylabel('BOLD (arb.units)');
    subplot(1,3,3) ; 
    hz = (freqs>80 & freqs<100) ;
    x = squeeze(mean(mean(mmtersp(:,stims,hz),3),1)); 
    y = squeeze(mean(meanf(:,stims),1)); 
    P = polyfit(x,y,1);yfit = P(1)*x+P(2);   
    hz = find(freqs>80 & freqs<100 ) ;  plot(x,yfit,'r','LineWidth',2);hold on ; 
    errorbarxy(x,y,squeeze(std(mean(mmtersp(:,stims,hz),3),0,1))./sqrt(22),squeeze(std(meanf(:,stims),0,1)./sqrt(22)),{'o','k','k'}) ; 
    xlabel('BBG (db)') ; ylabel('BOLD (arb.units)') ;
    
    subplot(1,3,3) ; 
    hz = (freqs>80 & freqs<100) ;
    x = squeeze(mean(mbolds(:,stims),1)) ; 
    y = squeeze(mean(meanf(:,stims),1)); 
    P = polyfit(x,y,1);yfit = P(1)*x+P(2);   
    hz = find(freqs>80 & freqs<100 ) ;  plot(x,yfit,'r','LineWidth',2);hold on ; 
    errorbarxy(x,y,squeeze(std(mbolds(:,stims),0,1))./sqrt(22),squeeze(std(meanf(:,stims),0,1)./sqrt(22)),{'o','k','k'}) ; 
    xlabel('BBG (db)') ; ylabel('BOLD (arb.units)') ;
    
    
    
    
    
    
    clear meanf 
    for i=1:22 ; figure
        for j=1:6  
            a = squeeze(zscore(mean(allstimepochs(i,j,:,5:6),4),0,3)) ; 
            goods = find(abs(a)<2) ; 
            meanf(i,j,:) = squeeze(mean(allstimepochs(i,j,goods,:),3)) ; 
            subplot(2,3,j); bar(a) ; 
        end  
    end
    
    
    
    
    
    
    %%%% stats 
    mersp = squeeze(mean(mean(allersp(:,:,1:3,:,:,:),3),4)) ; 
    mtersp = squeeze(mean(mersp(:,:,:,times>0 & times<2),4)) ; 
    nbg = find(freqs>60 & freqs<80) ; bbg = find((freqs > 40 &freqs<55 )| freqs>80) ; 
    abeta = find(freqs>10 & freqs<25) ; 
    p = anova1(squeeze(mean(mtersp(:,[1,6],nbg),3)),[],'off')

    
    anova1(squeeze(mean(mtersp(:,[1,5],nbg),3)),[],'off')
    
    
    vec = [squeeze(mean(mean(mtersp(:,stims,abeta),1),3));squeeze(mean(mean(meanf(:,stims,5:6),1),3))] ; 
    vec = [squeeze(mean(mbolds(:,stims)));squeeze(mean(mean(meanf(:,stims,5:6),1),3))] ; 

    
   
    
%%% ersp    
for st=1:length(stims) ; 
gfreqs = find(freqs>8) ; trimfreqs = freqs(gfreqs) ; 
fbands = [20,40,60,80,100] ; freqinds = zeros(size(fbands)) ; 
for i=1:length(fbands) ; diffs = abs(freqs(gfreqs)-fbands(i)) ; freqinds(i) = find(diffs==min(diffs)) ; end
tbands = [0,1,2] ; tinds = [44,101,159] ; 
subplot(2,3,st) ; mcersp = squeeze(mean(mtersp(:,1,:,:),1)) ;  
imagesc(mcersp(gfreqs,:),[-2,2]) ; 
set(gca,'YTick',freqinds,'YTickLabel',fbands,'XTick',tinds,'XTickLabel',tbands) ; axis xy
xlabel('time(s)') ; ylabel('frequency(hz)') ; vline([tinds(1),tinds(3)],'k') ;  
end


%%% anovas
mttersp = squeeze(mean(mtersp(:,:,:,times>0 & times<2),4)) ; 
for i=1:60
   [p,anovatab,stats] = anova1(squeeze(mttersp(:,[1,2],i)),[],'off') ;  
   cfvals(i) = anovatab{2,5} ; 
   cpvals(i) = p ; 
   [p,anovatab,stats] = anova1(squeeze(mttersp(:,[1,6],i)),[],'off') ;  
   rfvals(i) = anovatab{2,5} ; 
   rpvals(i) = p ; 
end

plot(cfvals) ; hold on ; plot(rfvals,'r') ; 

%%% t-test
mttersp = squeeze(mean(mtersp(:,:,:,times>0 & times<2),4)) ; 
for i=1:60
   [H,p,CI,stats] = ttest(squeeze(mttersp(:,[1],i))',squeeze(mttersp(:,[2],i))') ;  
   ctvals(i) = stats.tstat ; 
   cpvals(i) = p ; 
   [H,p,CI,stats] = ttest(squeeze(mttersp(:,[1],i))',squeeze(mttersp(:,[6],i))') ;  
   rtvals(i) = stats.tstat ; 
   rpvals(i) = p ; 
end
plot(ctvals) ; hold on ; plot(rtvals,'r') ; 
sig_contrast = cpvals < 0.01 ; 
sig_rnd = rpvals < 0.01 ; 




%%% plot the significance shite
gfreqs = find(freqs>8) ; trimfreqs = freqs(gfreqs) ; 
fbands = [20,40,60,80,100] ; freqinds = zeros(size(fbands)) ; 
for i=1:length(fbands) ; diffs = abs(freqs(gfreqs)-fbands(i)) ; freqinds(i) = find(diffs==min(diffs)) ; end

f = figure; set(f,'Position',[100,100,1000,800]) ; 
% contrast
stims1 = [1,3,2] ; colors1 = {[1,0,0],[.6,.6,0],[0,.6,0]} ; 
gfreqs = find(freqs>8) ; subplot(2,2,1) ; 
for s=1:length(stims1) ; 
errorbar(squeeze(mean(mttersp(:,stims1(s),gfreqs),1))',squeeze(std(mttersp(:,stims1(s),gfreqs),0,1))'./sqrt(22),'Color',colors1{s},'LineWidth',1) ; hold on  ; 
xlim([0,length(gfreqs)+1]) ; 
text(find(sig_contrast(gfreqs)==1),ones(1,length(find(sig_contrast(gfreqs)==1)))*1.3,'*')  ; ylim([-2,1.5]) ; hline(0,'k') ; 
end
xlabel('frequency(hz)') ; ylabel('power(db)') ; set(gca,'XTick',freqinds,'XTickLabel',fbands) ; legend({'unperturbed','33%MC','5%MC'}) ; 
title('* p[mean(unperturbed) - mean(5%MC) != 0] < 0.01') ; 

% randomization
stims1 = [1,5,6] ; colors1 = {[1,0,0],[0,.6,.6],[0,0,.6]} ; 
gfreqs = find(freqs>8) ; subplot(2,2,2) ; 
for s=1:length(stims1) ; 
errorbar(squeeze(mean(mttersp(:,stims1(s),gfreqs),1))',squeeze(std(mttersp(:,stims1(s),gfreqs),0,1))'./sqrt(22),'Color',colors1{s},'LineWidth',1) ; hold on  ; 
xlim([0,length(gfreqs)+1]) ; 
text(find(sig_rnd(gfreqs)==1),ones(1,length(find(sig_rnd(gfreqs)==1)))*1.3,'*')  ; ylim([-2,1.5]) ; hline(0,'k') ; 
end
xlabel('frequency(hz)') ; ylabel('power(db)') ; set(gca,'XTick',freqinds,'XTickLabel',fbands) ;  legend({'unperturbed','10%SR','60%SR'}) ; 
title('* p[mean(unperturbed) - mean(60%SR) != 0] < 0.01') ; 

%%%%%%% FMRI shite
meanf = squeeze(mean(allstimepochs,3)) ; 
for i=1:10  ;
   [H,p,CI,stats] = ttest(squeeze(meanf(:,1,i)),squeeze(meanf(:,2,i))) ;  
   fctvals(i) = stats.tstat ; 
   fcpvals(i) = p ; 
   [H,p,CI,stats] = ttest(squeeze(meanf(:,1,i)),squeeze(meanf(:,6,i))) ;  
   frtvals(i) = stats.tstat ; 
   frpvals(i) = p ;     
end
sig_fcontrast = fcpvals < 0.01 ; 
sig_frnd = frpvals < 0.01 ;


% contrast
f = figure; set(f,'Position',[100,100,1000,800]) ; 
stims1 = [1,3,2] ; colors1 = {[1,0,0],[.6,.6,0],[0,.6,0]} ; 
subplot(2,2,1) ; mts = 2:10 ; 
for s=1:length(stims1) ; 
    errorbar(squeeze(mean(meanf(:,stims1(s),mts),1))',squeeze(std(meanf(:,stims1(s),mts),0,1))'./sqrt(22),'Color',colors1{s},'LineWidth',2) ; hold on  ; 
    text(find(sig_fcontrast(mts)==1),ones(1,length(find(sig_fcontrast(mts)==1)))*1.3,'*')  ; ylim([-.75,1.5]) ; 
    xlim([.5,9.5]) ; hline(0,'k') ; 
end
legend({'unperturbed','33%MC','5%MC'}) ; xlabel('times(s)') ; ylabel('task-rest(a.u)') ; set(gca,'XTick',1:length(mts),'XTickLabel',0:2:16) ;
title('* p[mean(unperturbed) - mean(5%MC) != 0] < 0.01') ; 

% randomization
stims1 = [1,5,6] ; colors1 = {[1,0,0],[0,.6,.6],[0,0,.6]} ; 
subplot(2,2,2) ; mts = 2:10 ; 
for s=1:length(stims1) ; 
    errorbar(squeeze(mean(meanf(:,stims1(s),mts),1))',squeeze(std(meanf(:,stims1(s),mts),0,1))'./sqrt(22),'Color',colors1{s},'LineWidth',2) ; hold on  ; 
    text(find(sig_frnd(mts)==1),ones(1,length(find(sig_frnd(mts)==1)))*1.6,'*')  ; ylim([-.75,1.8]) ; 
    xlim([.5,9.5]) ; hline(0,'k') ; 
end
legend({'unperturbed','10%SR','60%SR'}) ; xlabel('times(s)') ; ylabel('task-rest(a.u)') ; set(gca,'XTick',1:length(mts),'XTickLabel',0:2:16) ;
title('* p[mean(unperturbed) - mean(60%SR) != 0] < 0.01') ; 


barwitherr(squeeze(std(mean(mttersp(:,[1,5,6],freqs>60 & freqs<80),3),0,1))./sqrt(22),squeeze(mean(mean(mttersp(:,[1,5,6],freqs>60 & freqs<80),1),3))) ;







%%%% single trial correlations
plot(squeeze(mean(mean(me(1,:,freqs>60 & freqs<80,times>0  & times<2),3),4)),squeeze(mean(mean(me(1,:,freqs>60 & freqs<90,times>0  & times<2),3),4)),'o') ; 

for i=1:60
    for j=1:60
        corrts(i,j) = corr2(squeeze(mean(me(1,i,times>0 & times<2),4)),squeeze(mean(me(1,i,times>0 & times<2),4))) ; 
        
    end
end

% bar charts
colors = [[0,.6,0];[.6,.6,0];[1,0,0];[0,.6,.6];[0,0,.6]] ;
%%%
mttersp = squeeze(mean(mtersp(:,:,:,times>0&times<2),4)) ; 
meanf = squeeze(mean(allstimepochs,3)) ; 
otherstims = [2,3,5,6] ; 

%%% alpha/beta
fHand = figure; set(fHand,'Position',[500,500,600,400]) ; 
y = squeeze(mean(mean(mttersp(:,stims,freqs>10 & freqs<25),3),1)) ;
s = squeeze(std(mean(mttersp(:,stims,freqs>10 & freqs<25),3),0,1))./sqrt(22) ;
aHand = axes('parent', fHand);hold(aHand, 'on')
for i = 1:numel(y);bar(i, y(i), 'parent', aHand, 'facecolor', colors(i,:));end
set(gca, 'XTick', 1:numel(y), 'XTickLabel', {'R0', 'R1', 'R2'})
h1 = errorbar(y,s,'k');set(h1,'linestyle','none')
ylabel('EEG 10-25Hz power(db)') ; ylim([-1.6,0]) ; 
for i=1:length(otherstims);[h(i),p(i),ci,stats] = ttest(squeeze(mean(mttersp(:,otherstims(i),freqs>10 & freqs<25),3))',squeeze(mean(mttersp(:,1,freqs>10 & freqs<25),3))') ; end
for i=1:length(p) ; ps{i} = ['p=',num2str(p(i))] ; end
ys = double(y([1,2,4,5])-.2) ; 
%ys=[.5,.5,.5,.5];
h=text([1,2,4,5]-.35,ys,ps);
title('alpha/beta p-values relative to unperturbed') ; 

%%% NBG
fHand = figure; set(fHand,'Position',[500,500,600,400]) ;  
y = squeeze(mean(mean(mttersp(:,stims,freqs>60 & freqs<80),3),1)) ;
s = squeeze(std(mean(mttersp(:,stims,freqs>60 & freqs<80),3),0,1))./sqrt(22) ;
aHand = axes('parent', fHand);hold(aHand, 'on')
for i = 1:numel(y);bar(i, y(i), 'parent', aHand, 'facecolor', colors(i,:));end
set(gca, 'XTick', 1:numel(y), 'XTickLabel', {'R0', 'R1', 'R2'})
h1 = errorbar(y,s,'k');set(h1,'linestyle','none')
ylabel('EEG 60-80Hz power(db)') ; ylim([0,1]) ; 
for i=1:length(otherstims);[h(i),p(i),ci,stats] = ttest(squeeze(mean(mttersp(:,otherstims(i),freqs>60 & freqs<80),3)),squeeze(mean(mttersp(:,1,freqs>60 & freqs<80),3))) ; end
for i=1:length(p) ; ps{i} = ['p=',num2str(p(i))] ; end
ys = double(y([1,2,4,5])+.13) ; 
%ys=[.5,.5,.5,.5];
h=text([1,2,4,5]-.35,ys,ps);
title('NBG p-values relative to unperturbed') ; 

%%% classical gamma
fHand = figure; set(fHand,'Position',[500,500,600,400]) ; 
y = squeeze(mean(mean(mttersp(:,stims,freqs>40 & freqs<90),3),1)) ;
s = squeeze(std(mean(mttersp(:,stims,freqs>40 & freqs<90),3),0,1))./sqrt(22) ;
aHand = axes('parent', fHand);hold(aHand, 'on')
for i = 1:numel(y);bar(i, y(i), 'parent', aHand, 'facecolor', colors(i,:));end
set(gca, 'XTick', 1:numel(y), 'XTickLabel', {'R0', 'R1', 'R2'})
h1 = errorbar(y,s,'k');set(h1,'linestyle','none')
ylabel('EEG 40-90Hz power(db)') ; 
for i=1:length(otherstims);[h(i),p(i),ci,stats] = ttest(squeeze(mean(mttersp(:,otherstims(i),freqs>40 & freqs<90),3)),squeeze(mean(mttersp(:,1,freqs>40 & freqs<90),3))) ; end
for i=1:length(p) ; ps{i} = ['p=',num2str(p(i))] ; end
ys = double(y([1,2,4,5])+.13) ; 
%ys=[.5,.5,.5,.5];
h=text([1,2,4,5]-.35,ys,ps);
title('gamma p-values relative to unperturbed') ; 

%%% BOLD
fHand = figure; set(fHand,'Position',[500,500,600,400]) ; 
y = squeeze(mean(mean(meanf(:,stims,5:6),1),3)) ;
s = squeeze(std(mean(meanf(:,stims,5:6),3),0,1))./sqrt(22) ;
aHand = axes('parent', fHand);hold(aHand, 'on')
for i = 1:numel(y);bar(i, y(i), 'parent', aHand, 'facecolor', colors(i,:));end
set(gca, 'XTick', 1:numel(y), 'XTickLabel', {'R0', 'R1', 'R2'})
h1 = errorbar(y,s,'k');set(h1,'linestyle','none') ; ylim([0,1.6]) ; 
ylabel('BOLD 6-8s task-rest') ; 
for i=1:length(otherstims);[h(i),p(i),ci,stats] = ttest(squeeze(mean(meanf(:,otherstims(i),5:6),3))',squeeze(mean(meanf(:,1,5:6),3))') ; end
for i=1:length(p) ; ps{i} = ['p=',num2str(p(i))] ; end
ys = double(y([1,2,4,5])+.13) ; 
%ys=[.5,.5,.5,.5];
h=text([1,2,4,5]-.35,ys,ps);
title('BOLD p-values relative to unperturbed') ; 







clear corrs
for s=1:22
meanfreqs = squeeze(mean(mean(allersp(s,:,1,:,:,times>0 & times<2),3),6)) ; 
for i=1:6
    for j=1:60
        for k=1:60
            corrs(s,i,j,k) = corr2(squeeze(meanfreqs(i,:,j)),squeeze(meanfreqs(i,:,k))) ; 
            
            
        end
    end
end
end
for i=1:22 ; figure ; for j=1:6 ; subplot(2,3,j) ; imagesc(squeeze(corrs(i,j,:,:)),[-1,1]) ; end ; end


responders = {'alexandra3','fabio','gab','gabriella','genevieve','gina','jeremie','julie','katrine','marie','maxime','mingham','po','russell','suhan2','tegan2'} ;


%%%%%%%%% compare with freesurfer labels
labs = parse_freesurfer ; 
meaneeg = squeeze(mean(mean(mtersp(:,:,:,times>0&times<2),2),4)) ; 
%meanf = squeeze(mean(mean(allstimepochs(:,:,:,5:6),3),4)) ; meaneeg = meanf ; 
clear corrmat ; 
for i=1:size(labs,2)
    labi = cell2mat({labs{:,i,2}}) ;
    for j=1:size(labi,1) ; 
        for k=1:size(meaneeg,2)
            corrmat{i}(j,k) = corr2(meaneeg(:,k),labi(j,:)') ;
        end
    end
end
for i=1:19 ; subplot(4,5,i) ; imagesc(squeeze(corrmat{i}),[-1,1]) ; end


imagesc(corrmat{9},[-1,1]) ; set(gca,'YTick',1:size(corrmat{1},1),'YTickLabel',labs{1,9,1}) ; 

param =9; hzind = 5; paramrow = 13 ; 
labi = cell2mat({labs{:,param,2}}) ;
figure,plot(squeeze(meaneeg(:,hzind)),labi(paramrow,:),'o') ; title(num2str(corr2(squeeze(meaneeg(:,hzind))',labi(paramrow,:)))) ; 
figure,imagesc(corrmat{param},[-1,1]) ; set(gca,'YTick',1:size(corrmat{param},1),'YTickLabel',labs{1,param,1}) ; 












