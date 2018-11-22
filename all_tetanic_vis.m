clear all ; close all ; 
stims = {'S  1','S  2','S 11','S 12'} ;
subs = {'AD','AM','CV','DT','FMP','KC','LB','RB','TD'};
bades = {[],[],[21,61],[9,42],[19,57,62],[1,32,33,37],[10],[],[18,50]};
subcomps = {[3,7],[3,4,6,19],[7,10,13,24],[4,10],[3,9,26],[4,5],[32],[3,17],[5,12]};

%{
labs1 = {'FP1','FZ','F3','F7','F5','FC5','FC1','C3','T7','TP9','CP5','CP1','PZ','P3','P7','O1','OZ','O2','P4','P8','TP10','CP6','CP2','CZ','C4','T8','FT10','FC6','FC2','F4','F8','FP2',...
    'PO8','AF7','AF3','AFZ','F1','F5','FT7','FC3','C1','C5','TP7','CP3','P1','P5','PO7','PO3','POZ','PO4','PO8','P6','P2','CPZ','CP4','TP8','C6','C2','FC4','FT8','F6','AF8','AF4','FZ'};

labs2 = {'FP1','FP2','F7','F3','FZ','F4','F8','FC5','FC1','FC2','FC6','T7','C3','CZ','C4','T8','TP9','CP5','CP1','CP2','CP6','TP10','P7','P3','PZ','P4','P8','PO9','O1','OZ','O2','PO10',...
    'IZ','AF7','AF3','AF4','AF8','F5','F1','F2','F6','FT9','FT7','FC3','FC4','FT8','FT10','C5','C1','C2','C6','TP7','CP3','CPZ','CP4','TP8','P5','P1','P2','P6','PO7','PO3','POZ','PO4'};
%}

cap = {'fp1','fz','f3','f7','ft9','fc5','fc1','c3','t7','tp9','cp5','cp1','pz','p3','p7','o1','oz','o2','p4','p8','tp10','cp6','cp2','cz','c4','t8','ft10','fc6','fc2','f4','f8','fp2',...
       'af7','af3','afz','f1','f5','ft7','fc3','c1','c5','tp7','cp3','p1','p5','po7','po3','poz','po4','po8','p6','p2','cpz','cp4','tp8','c6','c2','fc4','ft8','f6','af8','af4','f2','iz'};

jly = {'fp1','fp2','f3','f4','c3','c4','p3','p4','o1','o2','f7','f8','t7','t8','p7','p8','fz','cz','pz','oz','fc1','fc2','cp1','cp2','fc5','fc6','cp5','cp6','tp9','tp10','poz','ecg',...
    'f1','f2','c1','c2','p1','p2','af3','af4','fc3','fc4','cp3','cp4','po3','po4','f5','f6','c5','c6','p5','p6','af7','af8','ft7','ft8','tp7','tp8','po7','po8','ft9','ft10','fpz','cpz'};
   
final = {'fp1','fp2','f7','f3','fz','f4','f8','fc5','fc1','fc2','fc6','t7','c3','cz','c4','t8','tp9','cp5','cp1','cp2','cp6','tp10','p7','p3','pz','p4','p8','po9','o1','oz','o2','po10',...
    'af7','af3','af4','af8','f5','f1','f2','f6','ft9','ft7','fc3','fc4','ft8','ft10','c5','c1','c2','c6','tp7','cp3','cpz','cp4','tp8','p5','p1','p2','p6','po7','po3','poz','po4','po8'};
   

% need to move the final indices to the cap indices - find the 

cd(['E:\jf_data\Tetanic Visual\RB']) ; rmerged = pop_loadset('merged.set'); 

for sb=1:length(subs)
    
    cd(['E:\jf_data\Tetanic Visual\',subs{sb}]) ;
    
    %{
    discr = dir('*vhdr') ;
    for i=1:length(discr)
       EEG = pop_loadbv('.',discr(i).name) ; 
       EEG = pop_resample(EEG,256) ; 
       if i==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG) ; end
    end
    merged = pop_chanedit(merged,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
    merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61) - eegfiltfft(merged.data,merged.srate,84,86) - eegfiltfft(merged.data,merged.srate,0,1);  

    figure,bar(sum(abs(diff(merged.data,1,2)),2))
    %merged = pop_interp(merged,[32],'spherical');
    pop_saveset(merged,'merged.set'); 

    %}
    merged = pop_loadset('merged.set');
    %figure,bar(sum(abs(diff(merged.data,1,2)),2)) ; title(subs{sb}); 
   
    for i=1:length(bades{sb})
       merged.data(bades{sb}(i),:) = rand(1,length(merged.data))/100;  
        
    end
    
      %{
    mergefilt = merged ; mergefilt.data = eegfiltfft(merged.data,merged.srate,1,90)  ; 
    [weights,sphere] = runica(mergefilt.data(:,1:end),'maxsteps',128); 
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
    
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(allersp([1,3],i,:,:))),[-8,8]) ; colormap jet ;  title(i) ; axis xy ; end ;suptitle(subs{sb}); 
    figure,for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),rmerged.chanlocs) ; title(i); end ; suptitle(subs{sb}); 

    %}
    %{
    icaw = load('icaw'); icaw= icaw.icaw; 
    weights = icaw{1} ; sphere = icaw{2}; 
    winv = pinv(weights*sphere); 
    
    cap = {
        'fp1','fz','f3','f7','ft9','fc5','fc1','c3','t7','tp9','cp5','cp1','pz','p3','p7','o1','oz','o2','p4','p8','tp10','cp6','cp2','cz','c4','t8','ft10','fc6','fc2','f4','f8','fp2',...
       'af7','af3','afz','f1','f5','ft7','fc3','c1','c5','tp7','cp3','p1','p5','po7','po3','poz','po4','po8','p6','p2','cpz','cp4','tp8','c6','c2','fc4','ft8','f6','af8','af4','f2','iz',...
        };

    final = { 'fp1','fp2','f7','f3','fz','f4','f8','fc5','fc1','fc2','fc6','t7','c3','cz','c4','t8','tp9','cp5','cp1','cp2','cp6','tp10','p7','p3','pz','p4','p8','po9','o1','oz','o2','po10',...
        'af7','af3','af4','af8','f5','f1','f2','f6','ft9','ft7','fc3','fc4','ft8','ft10','c5','c1','c2','c6','tp7','cp3','cpz','cp4','tp8','p5','p1','p2','p6','po7','po3','poz','po4','po8',...
       
        };

    new_winv = winv;
    
    neweeg = merged; 
    %}
  
    %{
    for i=1:length(final)
        final_ind = find(strcmpi(final{i},cap)); 
        if ~isempty(final_ind)
             new_winv(:,i) = winv(:,final_ind); 
            %final_ind = find(strcmpi(final{final_ind},final)); 

            neweeg.data(i,:) = merged.data(final_ind,:); 
        else 
            neweeg.data(i,:) = 0; 
        end
    end
 
    %temp = new_winv(1:32,:); new_winv(1:32,:) = new_winv(33:64,:); new_winv(33:64,:) = temp; 
    figure,for i=1:16 ; subplot(5,5,i) ; topoplot(new_winv(:,i),merged.chanlocs) ; title(i); end
    %}
    
    
    %{
    [cweights,csphere] = runica(merged.data,'maxsteps',128); 
    cwinv = pinv(cweights*csphere); 
        figure,for i=1:16 ; subplot(5,5,i) ; topoplot(cwinv(:,i),merged.chanlocs) ; title(i); end
    
    %}
     icaw = load('icaw'); icaw= icaw.icaw; 
    weights = icaw{1} ; sphere = icaw{2}; 
    winv = pinv(weights*sphere); 
    comps = subcomps{sb}; 
    newmerged = merged ; newmerged.data = weights*sphere*merged.data ; 
    clear allersp erps
    for s=1:length(stims) ; disp(s) ; 
        allep = pop_epoch(newmerged,{stims{s}},[-1,3]) ; 
    for i=1:length(comps) 
        for t=1:size(allep.data,3)
            [allersp(s,i,t,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(comps(i),:,t)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',NaN,'verbose','off','timesout',100) ; 
        end
        erps(s,i,:,1:size(allep.data,3)) = allep.data(comps(i),:,1:end); 
    end
    end

    bersp = allersp - repmat(mean(allersp(:,:,:,:,times<0),5),[1,1,1,1,100]); 
    sb_bersp(sb,:,:,:,:) = squeeze(mean(mean(bersp,2),3)); 

    
    %{
    figure,
    trigtitles = {'S  1 (pre @ tetanic freq)','S  2 (pre @ lower freq)','S 11 (post @ tetanic freq)','S 12 (post @ lower freq)'};
    for i=1:4 
        subplot(3,2,i);imagesc(times,freqs,squeeze(mean(mean(bersp(i,1,:,:,:),2),3)),[-7,7]); axis xy ; colormap jet; xlabel('time(ms)') ; ylabel('frequency(hz)'); title(trigtitles{i}); vline([0,2],'k'); vline(0.5,'m');
    end
    
    sb_bersp(sb,:,:,:,:) = squeeze(mean(mean(bersp,2),3)); 
    %}
end


titles = {'PRE high freq','PRE low freq','POST high freq','POST low freq'};
figure,
for i=1:4 ; subplot(2,5,i); imagesc(times,freqs,squeeze(mean(sb_bersp(:,i,:,:),1)),[-3,3])  ; axis xy; colormap jet ; title(titles{i}); xlabel('time(s)') ; ylabel('frequency(hz)');  end
subplot(2,5,6); shadedErrorBar(freqs,squeeze(mean(mean(sb_bersp(:,[3],:,times>0 & times<2),1),4))',squeeze(std(mean(sb_bersp(:,[3],:,times>0 & times<2),4),0,1))'/sqrt(5),{'b'}) ; ylim([-5,3]); hold on; 
shadedErrorBar(freqs,squeeze(mean(mean(sb_bersp(:,[1],:,times>0 & times<2),1),4))',squeeze(std(mean(sb_bersp(:,[1],:,times>0 & times<2),4),0,1))'/sqrt(5),{'r'}) ; ylim([-5,3]);
title('high freq'); xlim([0,120]);
subplot(2,5,7); shadedErrorBar(freqs,squeeze(mean(mean(sb_bersp(:,[2],:,times>0 & times<2),1),4))',squeeze(std(mean(sb_bersp(:,[2],:,times>0 & times<2),4),0,1))'/sqrt(5),{'b'}) ; ylim([-5,3]); hold on; 
shadedErrorBar(freqs,squeeze(mean(mean(sb_bersp(:,[4],:,times>0 & times<2),1),4))',squeeze(std(mean(sb_bersp(:,[4],:,times>0 & times<2),4),0,1))'/sqrt(5),{'r'}) ; ylim([-5,3]); xlabel('frequency(hz)'); ylabel('db'); 
title('low freq'); xlim([0,120]);

mbersp = squeeze(mean(mean(sb_bersp(:,:,25:30,times>0 & times<2),3),4)); 
subplot(2,5,8);
barwitherr(squeeze(std(mbersp(:,[1,3]),0,1))/3,mean(mbersp(:,[1,3]),1)); ylabel('dB'); set(gca,'XTickLabel',{'pre','post'});
[h,p,ci,stats] = ttest(mbersp(:,3),mbersp(:,1)); title(['high freq, ',format_t(stats.tstat),' ',format_p(p)]); 
subplot(2,5,9); 
barwitherr(squeeze(std(mbersp(:,[2,4]),0,1))/3,mean(mbersp(:,[2,4]),1)); ylabel('dB'); set(gca,'XTickLabel',{'pre','post'}); 
[h,p,ci,stats] = ttest(mbersp(:,4),mbersp(:,2)); title(['low freq, ', format_t(stats.tstat),' ',format_p(p)]); 
subplot(2,5,10);
plot(1,'b') ; hold on ; plot(2,'r'); legend({'post','pre'});


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

%{
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
%}