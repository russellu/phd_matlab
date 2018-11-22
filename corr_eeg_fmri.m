clear all ; close all

subs = {'charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','lisa','marc','marie',...
    'mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ;


clear allersp ; 
for s=1:length(subs)
    cd(['c:/shared/allres/',subs{s}]) ; ls 

    ersp = load('ersp.mat') ; 
    stersp = load('stersp') ; stersp = stersp.ersp ; 
    allstersp(s,:,:,:,:) = squeeze(stersp(:,1:3,:,:)) ;
    ersp = ersp.ersp ; allersp(s,:,:) = squeeze(mean(ersp(1:2,:,:),1)) ; 
    
    %subplot(4,6,s) ; imagesc(squeeze(mean(ersp(:,:,:))),[-2,2]) ; 
    cd(['C:\shared\all_white_normals\fmris\sub_',subs{s}]) ; 
    pchangei = load('pchangei') ; pchangei = pchangei.pchangei ; 
    allpchange(s) = pchangei ; 
    pvoxi = load('pvoxi') ; pvoxi = pvoxi.pvoxi ; 
    allpvox(s) = pvoxi ; 
    posi = load('allposi') ; posi = posi.allposi ; 
    allposi(s,:,:,:,:,:) = posi ; 
    
end

fpercs = squeeze(mean(mean(allposi,2),4)) ; 
fpercs = (squeeze(mean(fpercs(:,:,:,5:6),4)) - squeeze(mean(fpercs(:,:,:,1:2),4))) ./ squeeze(mean(fpercs(:,:,:,1:2),4)) ; 
fpercs = squeeze(mean(fpercs,3)) ; 
fr = corr(fpercs) ; 

mersp = squeeze(mean(mean(allstersp(:,:,:,:,40:180),3),5)) ; 

for i=1:6
    for j=1:60
        corrs(i,j) = corr2(squeeze(mersp(:,i,j)),fpercs(:,i)) ; 
    end
end
for i=1:60 
   subcorrs(i) = corr2(squeeze(mean(mersp(:,:,i),2)),squeeze(mean(fpercs,2))) ;  
end


for i=1:22
    for j=1:60
        scorrs(i,j) = corr2(squeeze(mersp(i,:,j)),fpercs(i,:)) ; 
    end
end
for i=1:60 
   stimcorrs(i) = corr2(squeeze(mean(mersp(:,:,i),1)),squeeze(mean(fpercs,1))) ;  
end

for i=1:6
    for j=1:60
        voxcorrs(i,j) = corr2(allpvox,squeeze(mersp(:,i,j))') ; 
    end
end
for i=1:60
    meanvoxcorrs(i) = corr2(allpvox,squeeze(mean(mersp(:,:,i),2))') ; 
end




subplot(2,2,3) ; 
shadedErrorBar([],mean(voxcorrs,1),std(voxcorrs,0,1)/sqrt(6)) ; hold on ; plot(meanvoxcorrs,'r','LineWidth',2) ; 
set(gca,'XTick',1:5:60,'XTickLabel',(1:5:60)*2) ; xlabel('frequency(hz)') ; ylabel('correlation (pearson r)') ; 
title('subject specific correlations (EEG modulation vs #active voxels') ; hline(0,'k') ; 

subplot(2,2,1) ; 
shadedErrorBar([],mean(corrs,1),std(corrs,0,1)/sqrt(6)) ; hline(0,'k') ; hold on ; plot(subcorrs,'r','LineWidth',2) ; 
set(gca,'XTick',1:5:60,'XTickLabel',(1:5:60)*2) ; xlabel('frequency(hz)') ; ylabel('correlation (pearson r)') ; 
title('subject specific correlations (EEG modulation vs BOLD %change)') ; 
subplot(2,2,2) ; 
shadedErrorBar([],mean(scorrs,1),std(scorrs,0,1)/sqrt(22)) ; hline(0,'k') ; hold on ; plot(stimcorrs,'r','LineWidth',2) ; 
set(gca,'XTick',1:5:60,'XTickLabel',(1:5:60)*2) ; xlabel('frequency(hz)') ; ylabel('correlation (pearson r)') ; 
title('stimulus specific correlations (EEG modulation vs BOLD %change)') ; 


plot(1,'r') ; hold on ; plot(2,'k') ; legend({'correlation of average','average of correlations'}) ; 






