clear all ; close all ; 
stims = {'S  1','S  2','S 11','S 12'} ;
bades = {[32]}; 
%{
labs1 = {'FP1','FZ','F3','F7','F5','FC5','FC1','C3','T7','TP9','CP5','CP1','PZ','P3','P7','O1','OZ','O2','P4','P8','TP10','CP6','CP2','CZ','C4','T8','FT10','FC6','FC2','F4','F8','FP2',...
    'PO8','AF7','AF3','AFZ','F1','F5','FT7','FC3','C1','C5','TP7','CP3','P1','P5','PO7','PO3','POZ','PO4','PO8','P6','P2','CPZ','CP4','TP8','C6','C2','FC4','FT8','F6','AF8','AF4','FZ'};

labs2 = {'FP1','FP2','F7','F3','FZ','F4','F8','FC5','FC1','FC2','FC6','T7','C3','CZ','C4','T8','TP9','CP5','CP1','CP2','CP6','TP10','P7','P3','PZ','P4','P8','PO9','O1','OZ','O2','PO10',...
    'IZ','AF7','AF3','AF4','AF8','F5','F1','F2','F6','FT9','FT7','FC3','FC4','FT8','FT10','C5','C1','C2','C6','TP7','CP3','CPZ','CP4','TP8','P5','P1','P2','P6','PO7','PO3','POZ','PO4'};
%}

for sub=1%:length(subs)
    
    cd(['E:\jf_data\Tetanic Visual']) ;  
    discr = dir('C1_FMP_Tetanic Visuel.vhdr') ;
    for i=1:length(discr)
       EEG = pop_loadbv('.',discr(i).name) ; 
       EEG = pop_resample(EEG,256) ; 
       if i==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG) ; end
    end
    merged = pop_chanedit(merged,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
    merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61) - eegfiltfft(merged.data,merged.srate,84,86) ;  
    
    figure,bar(sum(abs(diff(merged.data,1,2)),2))
    merged = pop_interp(merged,[32],'spherical');
   
    mergefilt = merged ; mergefilt.data = eegfiltfft(merged.data,merged.srate,1,90)  ; 
    [weights,sphere] = runica(mergefilt.data(:,1:3:end),'maxsteps',128); 
    ica = merged; 
    ica.icaact = weights*sphere*merged.data ; winv = pinv(weights*sphere) ; 
    icaw{1} = weights; icaw{2} = sphere ; save('icaw','icaw') ; 
    winv = pinv(weights*sphere); 
    ica.data = ica.icaact ; clear allersp ; 
    newmerged = merged ; newmerged.data = weights*sphere*merged.data ; 
    for s=1:length(stims) ; disp(s) ; 
        allep = pop_epoch(newmerged,{stims{s}},[-1,3]) ; 
    for i=1:64 
            [allersp(s,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',0,'verbose','off','timesout',200) ; 
    end
    end
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(allersp([1,3],i,:,:))),[-8,8]) ; colormap jet ;  title(i) ; axis xy ; end
   
end
%{
elabs = {EEG.chanlocs.labels};
cwinv = winv; 
elab_inds = zeros(1,64); 
for i=1:size(elab_inds,2)
    if ~isempty(find(strcmpi(labs1{i},elabs)))
    elab_inds(i) = find(strcmpi(labs1{i},elabs)); 
    end
    
end
for i=1:size(elab_inds)
   if elab_inds(i) ~=0 
      cwinv(elab_inds(i),:) = winv(i,:);  
   end   
end
%}

comps = [3,4,6]; 
clear allersp erps
for s=1:length(stims) ; disp(s) ; 
    allep = pop_epoch(newmerged,{stims{s}},[-1,3]) ; 
for i=1:length(comps) 
    for t=1:30
        [allersp(s,i,t,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(comps(i),:,t)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',NaN,'verbose','off','timesout',100) ; 
    end
    erps(s,i,:,:) = allep.data(comps(i),:,1:30); 
end
end

bersp = allersp - repmat(mean(allersp(:,:,:,:,times<0),5),[1,1,1,1,100]); 

trigtitles = {'S  1 (pre @ tetanic freq)','S  2 (pre @ lower freq)','S 11 (post @ tetanic freq)','S 12 (post @ lower freq)'};
for i=1:4 
    subplot(3,2,i);imagesc(times,freqs,squeeze(mean(mean(bersp(i,1,:,:,:),2),3)),[-7,7]); axis xy ; colormap jet; xlabel('time(ms)') ; ylabel('frequency(hz)'); title(trigtitles{i}); vline([0,2],'k'); vline(0.5,'m');
end

for i=1:60 ; [h,p,ci,stats] = ttest2(squeeze(mean(bersp(3,1,:,i,times>0.5 & times<2),5)),squeeze(mean(bersp(1,1,:,i,times>0.5 & times<2),5))); ts1(i) = stats.tstat; ps1(i) = p;  end
for i=1:60 ; [h,p,ci,stats] = ttest2(squeeze(mean(bersp(4,1,:,i,times>0.5 & times<2),5)),squeeze(mean(bersp(2,1,:,i,times>0.5 & times<2),5))); ts2(i) = stats.tstat; ps2(i) = p;  end

subplot(3,2,5); 
shadedErrorBar(freqs,squeeze(mean(mean(bersp(1,1,:,:,times>0.5 & times<2),3),5)),squeeze(std(mean(bersp(1,1,:,:,times>0.5 & times<2),5),0,3))/sqrt(50),{'b'}); hold on ; xlim([0,120]);
shadedErrorBar(freqs,squeeze(mean(mean(bersp(3,1,:,:,times>0.5 & times<2),3),5)),squeeze(std(mean(bersp(2,1,:,:,times>0.5 & times<2),5),0,3))/sqrt(50),{'r'}); xlabel('frequency(hz)'); ylabel('db'); hline(0,'k'); ylim([-6,8]);
for i=1:length(ps1) ; if ps1(i) < 0.05 ; text(freqs(i),7,'*'); end ;end 
title('tetanic freq')

subplot(3,2,6); 
shadedErrorBar(freqs,squeeze(mean(mean(bersp(2,1,:,:,times>0.5 & times<2),3),5)),squeeze(std(mean(bersp(1,1,:,:,times>0.5 & times<2),5),0,3))/sqrt(50),{'b'}); hold on ; xlim([0,120]);
shadedErrorBar(freqs,squeeze(mean(mean(bersp(4,1,:,:,times>0.5 & times<2),3),5)),squeeze(std(mean(bersp(2,1,:,:,times>0.5 & times<2),5),0,3))/sqrt(50),{'r'}); xlabel('frequency(hz)'); ylabel('db'); hline(0,'k');ylim([-6,8]); 
for i=1:length(ps1) ; if ps2(i) < 0.05 ; text(freqs(i),7,'*'); end ;end
title('low freq');


% erp analysis
clear allersp
for s=1:length(stims) ; disp(s) ; 
    allep = pop_epoch(mergefilt,{stims{s}},[-1,3]) ; 
    for i=1:64    
        erps(s,i,:,:) = allep.data(i,:,:); 
    end
end
es = [61,62,63,64,29,30,31];

merps = squeeze(mean(erps(:,es,:,:),2)); 

for i=1:1024 ; [h,p,ci,stats] = ttest2(squeeze(merps(3,i,:)),squeeze(merps(1,i,:))) ; e_ps1(i) = p ; end
for i=1:1024 ; [h,p,ci,stats] = ttest2(squeeze(merps(4,i,:)),squeeze(merps(2,i,:))) ; e_ps2(i) = p ; end

subplot(2,1,1); 
shadedErrorBar(allep.times,squeeze(mean(mean(erps(1,es,:,:),4),2)),squeeze(std(mean(erps(1,es,:,:),2),0,4))/sqrt(50),{'b'}) ; hold on ; 
shadedErrorBar(allep.times,squeeze(mean(mean(erps(3,es,:,:),4),2)),squeeze(std(mean(erps(3,es,:,:),2),0,4))/sqrt(50),{'r'}) ;  ylim([-10,20]); vline(0,'k');
for i=1:length(e_ps1) ; if e_ps1(i) < 0.001 ; text(allep.times(i),19,'*'); end ; end
title('pre vs post ERP tetanic freq'); 

subplot(2,1,2); 
shadedErrorBar(allep.times,squeeze(mean(mean(erps(2,es,:,:),4),2)),squeeze(std(mean(erps(2,es,:,:),2),0,4))/sqrt(50),{'b'}) ; hold on ; 
shadedErrorBar(allep.times,squeeze(mean(mean(erps(4,es,:,:),4),2)),squeeze(std(mean(erps(4,es,:,:),2),0,4))/sqrt(50),{'r'}) ;  ylim([-10,20]); vline(0,'k');
for i=1:length(e_ps1) ; if e_ps2(i) < 0.001 ; text(allep.times(i),19,'*'); end ; end
title('pre vs post ERP lower freq'); 

figure,
plot([1],'b') ; hold on ; plot(2,'r') ; legend({'pre tetanic','post tetanic'});
