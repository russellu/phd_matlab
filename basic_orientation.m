clear all ; close all; 
subs = {'anne_marie','felix','florence','greg','jf','joel','lyes','mael','nic','russ','valerie'}; 
bades = {[],[],[],[48],[],[],[],[],[],[],[]}; 
subcomps = {[4,9],[4,7],[3,5,6],[2,5,7],[4,7],[3,8,9,10,11],[11,13,20],[6,8],[5,18],[11,16,27],[6,8,15]};
stims = {'S  1','S  2','S  3','S  4'};
swaps = {0,0,0,1,0,0,0,0,0,0,0};

for sb=1:length(subs)
    disp(subs{sb}); 
    cd(['E:\jly_orientation\',subs{sb}]);
    %{
    vhdrs = dir('*vhdr'); 
    for vhdr = 1:length(vhdrs)
        eeg = pop_loadbv('.',vhdrs(vhdr).name); 
        eeg = pop_chanedit(eeg,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
        eeg = pop_resample(eeg,250); 
        if vhdr==1 ; merged = eeg ; else merged = pop_mergeset(merged,eeg); end 
    end
    
    %figure,bar(zscore(sum(abs(diff(merged.data,1,2)),2))); suptitle(subs{sb}); 
    merged = pop_interp(merged,bades{sb},'spherical'); 
    pop_saveset(merged,'merged.set'); 
    %}
    
    %{
    merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59.5,60.5) - eegfiltfft(merged.data,merged.srate,84.5,85.5); 
    merged=  pop_saveset(merged,'merged.set'); 
    %}
    
    merged = pop_loadset('merged.set'); 

    filtmerged = merged ; filtmerged.data = eegfiltfft(merged.data, merged.srate, 1, 125); 
    
    %ep_all = pop_epoch(filtmerged,{'S  1','S  2','S  3','S  4'},[-4,24]); 
    
    %resall = reshape(ep_all.data,[64,numel(ep_all.data(1,:,:))]); 
    
    %for i=1:size(ep_all.data,3); resall(:,i:i+size(ep_all.data,2)-1) = ep_all.data(:,:,i) ; end
    %{
    [weights,sphere] = runica(filtmerged.data,'maxsteps',128); 
    saveica{1} = weights; saveica{2} = sphere; 
    save('saveica','saveica'); 
    
    winv = pinv(weights*sphere); 
    clear allersp ; 
    newmerged = merged ; newmerged.data = weights*sphere*merged.data ; 
    for s=1:length(stims) ; disp(s) ; 
        allep = pop_epoch(newmerged,{stims{s}},[-2,20]) ; 
        [sv,si] = sort( zscore(mean(std(allep.data(:,:,:),0,2),1)),'descend'); 
    for i=1:64 
            [allersp(s,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,5:end)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',0,'verbose','off','timesout',200) ; 
    end
    end   
    
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(allersp([1,2,3,4],i,:,:))),[-4,4]) ; colormap jet ;  title(i) ; axis xy ; end ;suptitle(subs{sb}); 
    figure,for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),merged.chanlocs) ; title(i); end ; suptitle(subs{sb}); 
  %}
    
    
    
    saveica = load('saveica'); saveica = saveica.saveica; weights = saveica{1}; sphere = saveica{2}; 
    
    winv = pinv(weights*sphere); 
    clear allersp ; 
    newmerged = merged ; newmerged.data = weights*sphere*merged.data ; 
    for s=1:length(stims) ; disp(s) ; 
        allep = pop_epoch(newmerged,{stims{s}},[-2,20]) ; 
    for i=1:length(subcomps{sb}) 
        for j=1:size(allep.data,3)
            [allersp(s,i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(subcomps{sb}(i),:,j)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',130,'baseline',NaN,'verbose','off','timesout',600) ; 
        end
    end
    end   
    bersp = allersp - repmat(squeeze(mean(allersp(:,:,:,:,times<0),5)),[1,1,1,1,600]); 

    %{
    sis = load('sis.mat'); sis = sis.sis;  
    clear cleanersp ; 
    newmerged = merged ; newmerged.data = weights*sphere*merged.data ; 
    for s=1:length(stims) ; disp(s) ; 
        allep = pop_epoch(newmerged,{stims{s}},[-2,20]) ; 
        
    for i=1:length(subcomps{sb}) 
            [cleanersp(s,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(subcomps{sb}(i),:,sis(s,5:end))),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',130,'baseline',0,'verbose','off','timesout',600) ; 
    end
    end   
    %}
    time_inds = find(times > 0 & times < 18); 
    
    gtime_inds = find(times>0.5 & times<18); 
    
    tbersp = bersp(:,:,:,:,time_inds); 

    
    revrots = load('E:\jly_orientation\revrots') ; revrots = revrots.revrots ; 
    fwdrots = load('E:\jly_orientation\fwdrots') ; fwdrots = fwdrots.fwdrots ; 
    
    horizontal = [355:359,1:5,175:185];
    vertical = [85:95,265:275];
    oblique = [40:50,130:140,220:230,310:320];
    cardinal = [1:10,85:95,175:185,265:275];
    quads = [315,225,135,45];
    bottom = 225:315; top = 45:135;
    right = [315:359,1:45]; left = [135:225];

    s1_angles = [315:-1:0,359:-1:270];
    s2_angles = [90:-1:0,359:-1:45];
    s3_angles = [225:1:359,0:1:270];
    s4_angles = [45:1:359,0:1:90];

    tbersp = bersp(:,:,:,:,time_inds); 
    clear res_angles
    res_angles(1,:) = imresize(s1_angles,[1,length(tbersp)],'nearest'); 
    res_angles(2,:) = imresize(s2_angles,[1,length(tbersp)],'nearest'); 
    res_angles(3,:) = imresize(s3_angles,[1,length(tbersp)],'nearest'); 
    res_angles(4,:) = imresize(s4_angles,[1,length(tbersp)],'nearest'); 
    
    tbersp = bersp(:,:,:,:,gtime_inds); 
    res_angles = res_angles(:,(length(time_inds)-length(gtime_inds)+1):end); 
    
    
    
   % plot(squeeze(mean(mean(mean(compdb(:,[2,4],14:15,:),1),2),3))) ; vline(oblique,'r'); vline(cardinal,'g');

    tbersp_amp = zeros(size(tbersp,1),size(tbersp,2),size(tbersp,3),size(tbersp,4),359); 
    
    for i=0:359
        for j=1:4
            bersp_ind = find(i==res_angles(j,:)); 
            tbersp_amp(j,:,:,:,i+1) = squeeze(mean(tbersp(j,:,:,:,bersp_ind),5)); 
        end
    end
    save('tbersp_amp','tbersp_amp'); 
    clear sis; 
    for i=1:4
        [sv,si] = sort(squeeze(abs(zscore(mean(mean(mean(abs(tbersp_amp(i,:,:,:,:)),4),5),2),[],3))),'descend'); 
        sis(i,:) = si; 
        sub_db(sb,i,:,:) = squeeze(mean(mean(tbersp_amp(i,:,si(5:end),:,:),3),2)); 
    end
    save('sis','sis'); 
    
end
%for i=1:11 ; figure ; for j=1:64 ; subplot(5,13,j) ;
%imagesc(squeeze(mean(sub_db(i,:,j,:,:),2)),[-3,3]) ; colormap jet; end ;
%suptitle(subs{i}); end

compdb = sub_db; 

cd E:\saved
stims_wedge = load('stims_wedge'); stims_wedge = stims_wedge.stims_wedge; 
stims_bars = load('stims_bars'); stims_bars = stims_bars.stims_bars; 
wedge_labels = {'0\circ','45\circ','90\circ','135\circ','180\circ','225\circ','270\circ','315\circ'};
band_names = {'mid gamma (40-70Hz','low gamma (28-32Hz)','beta (15-25Hz)','alpha (8-12Hz)'}; 


% FIGURE 1
figure,
for i=1:8 
   subplottight(8,10,(i-1)*10+1) ; imshow((stims_wedge(:,:,i))) ; 
   text(310,150,wedge_labels(i)); 
end

figure,
subplot(2,6,1); 
imagesc(1:360,freqs,squeeze(mean(mean(compdb(:,[1,3],:,:),1),2))) ; 
axis xy ; colormap jet; h=colorbar; title(h,'dB'); 
xlabel('wedge location (deg\circ)'); ylabel('frequency(hz)'); title('EEG retinotopic tuning');

subplot(2,6,2); 
plot(1:360,squeeze(mean(mean(mean(compdb(:,[1,3],20:35,:),1),2),3)),'Color',[0,0,0.5],'LineWidth',1); hold on; 
plot(1:360,squeeze(mean(mean(mean(compdb(:,[1,3],14:15,:),1),2),3)),'Color',[0,0.5,0],'LineWidth',1); 
plot(1:360,squeeze(mean(mean(mean(compdb(:,[1,3],8:12,:),1),2),3)),'Color',[0.5,0,0],'LineWidth',1); 
plot(1:360,squeeze(mean(mean(mean(compdb(:,[1,3],4:6,:),1),2),3)),'Color',[0.25,0.15,0.15],'LineWidth',1); hline(0,'k'); xlim([0,360]); 
title('frequency specific retinotopic tuning'); xlabel('wedge location (deg\circ)'); ylabel('dB'); 
grid on ;  set(gca,'XTick',[0,45,90,135,180,225,270,315],'XTickLabel',[0,45,90,135,180,225,270,315]); 
%for i=1:length(quads) ; vline(quads(i),'r'); text(quads(i),0.8,quadtitles{i}); end
legend(band_names); ylim([-1.5,0.7])

for i=1:60 
    [h,p,ci,stats] = ttest(squeeze(mean(mean(compdb(:,[1,3],i,top),2),4)),squeeze(mean(mean(compdb(:,[1,3],i,bottom),2),4)));     
    topbot_tvals(i) = stats.tstat; topbot_pvals(i) = p; 
    [h,p,ci,stats] = ttest(squeeze(mean(mean(compdb(:,[1,3],i,left),2),4)),squeeze(mean(mean(compdb(:,[1,3],i,right),2),4))); 
    leftright_tvals(i) = stats.tstat; leftright_pvals(i) = p; 
end
subplot(2,6,3); 
shadedErrorBar(freqs,squeeze(mean(mean(mean(compdb(:,[1,3],:,top),1),2),4)),squeeze(std(mean(mean(compdb(:,[1,3],:,top),2),4),0,1))/sqrt(11),'k'); hline(0,'k') ; hold on ; 
shadedErrorBar(freqs,squeeze(mean(mean(mean(compdb(:,[1,3],:,bottom),1),2),4)),squeeze(std(mean(mean(compdb(:,[1,3],:,bottom),2),4),0,1))/sqrt(11),{'m'}); hline(0,'k') ; 
ylabel('dB'); xlabel('frequency(hz)'); xlim([0,120]); title('upper vs lower, *p<0.01'); ylim([-2.5,1.5])
for i=1:length(topbot_tvals) ; if topbot_pvals(i) < 0.01; text(freqs(i),1.4,'*'); end; end

%{
subplot(2,6,4); plot(freqs,topbot_tvals); hold on ; plot(freqs,leftright_tvals); ylim([-6,6]);
ylabel('t-value'); xlabel('frequency(hz)'); xlim([0,120]); hline(0,'k'); hline([2.15,-2.15],'r'); text(122,2.2,'p<0.05'); text(122,-2.2,'p<0.05'); 
legend({'t (upper-lower)','t (left-right)'}); title('t-test: upper-lower, left-right')
%}
subplot(2,6,6) ; plot(1,'k') ; hold on ; plot(2,'m'); legend({'upper visual field','lower visual field'}); 

bands = {find(freqs>=40 & freqs<=70),find(freqs>=28 & freqs<=32),find(freqs>=15 & freqs<=25),find(freqs>=8 & freqs<=12)};

bandcolors = {[0,0,0.5],[0,0.5,0],[0.5,0,0],[0.25,0.15,0.15]};
subplots = [21,22,23,24]; 
for i=1:length(bands)
    subplot(3,10,subplots(i)); 
    b1 = bar(1,squeeze(mean(mean(mean(mean(compdb(:,[1,3],bands{i},bottom))))))); hold on; 
    b1.FaceColor = bandcolors{i}; % b1.LineStyle = ':'; b1.LineWidth = 3; b1.EdgeColor = 'r'; 
    b2 = bar(2,squeeze(mean(mean(mean(mean(compdb(:,[1,3],bands{i},top))))))); 
    b2.FaceColor = bandcolors{i}; 
    
    [h,p,ci,stats] = ttest(squeeze(mean(mean(mean(compdb(:,[1,3],bands{i},bottom),2),3),4)),squeeze(mean(mean(mean(compdb(:,[1,3],bands{i},top),2),3),4))); 
    
    mean_ebar = [squeeze(mean(mean(mean(mean(compdb(:,[1,3],bands{i},bottom)))))),squeeze(mean(mean(mean(mean(compdb(:,[1,3],bands{i},top))))))]; 
    std_ebar = [squeeze(std(mean(mean(mean(compdb(:,[1,3],bands{i},bottom),2),3),4),0,1))/sqrt(11),squeeze(std(mean(mean(mean(compdb(:,[1,3],bands{i},top),2),3),4),0,1))/sqrt(11)]; 
    errorbar([1,2],mean_ebar,std_ebar,'k.'); if i==1;  ylabel('dB'); end
    set(gca,'XTick',[1,2],'XTickLabel',{'upper','lower'});  xlim([0.5,2.5]); xlabel(band_names{i}); title([format_t(stats.tstat),', ' format_p(p)]);
end
set(gcf,'Position',[100 100 2000 450])



% FIGURE 2
figure,
for i=1:8 
   subplottight(8,10,(i-1)*10+1) ; imshow((stims_bars(:,:,i))) ; 
   text(310,150,wedge_labels(i)); 
end

figure,
subplot(2,6,1); 
imagesc(1:360,freqs,squeeze(mean(mean(compdb(:,[2,4],:,:),1),2))) ; 
axis xy ; colormap jet; h=colorbar; title(h,'dB'); 
xlabel('bar angle (deg\circ)'); ylabel('frequency(hz)'); title('EEG orientation tuning');

subplot(2,6,2); 
plot(1:360,squeeze(mean(mean(mean(compdb(:,[2,4],20:35,:),1),2),3)),'Color',[0,0,0.5],'LineWidth',1); hold on; 
plot(1:360,squeeze(mean(mean(mean(compdb(:,[2,4],14:15,:),1),2),3)),'Color',[0,0.5,0],'LineWidth',1); 
plot(1:360,squeeze(mean(mean(mean(compdb(:,[2,4],8:12,:),1),2),3)),'Color',[0.5,0,0],'LineWidth',1); 
plot(1:360,squeeze(mean(mean(mean(compdb(:,[2,4],4:6,:),1),2),3)),'Color',[0.25,0.15,0.15],'LineWidth',1); hline(0,'k'); xlim([0,360]); 
title('frequency specific orientation tuning'); xlabel('bar orientation (deg\circ)'); ylabel('dB'); 
grid on ; set(gca,'XTick',[0,45,90,135,180,225,270,315],'XTickLabel',[0,45,90,135,180,225,270,315]); 
%for i=1:length(quads) ; vline(quads(i),'r'); text(quads(i),0.8,quadtitles{i}); end
legend(band_names); 
for i=1:60 
    [h,p,ci,stats] = ttest(squeeze(mean(mean(compdb(:,[2,4],i,cardinal),2),4)),squeeze(mean(mean(compdb(:,[2,4],i,oblique),2),4)));     
    cardob_tvals(i) = stats.tstat; cardob_pvals(i) = p; 
end
subplot(2,6,3); 
shadedErrorBar(freqs,squeeze(mean(mean(mean(compdb(:,[2,4],:,cardinal),1),2),4)),squeeze(std(mean(mean(compdb(:,[2,4],:,cardinal),2),4),0,1))/sqrt(11),'k'); hline(0,'k') ; hold on ; 
shadedErrorBar(freqs,squeeze(mean(mean(mean(compdb(:,[2,4],:,oblique),1),2),4)),squeeze(std(mean(mean(compdb(:,[2,4],:,oblique),2),4),0,1))/sqrt(11),{'m'}); hline(0,'k') ; 
ylabel('dB'); xlabel('frequency(hz)'); xlim([0,120]); title('cardinal vs oblique, *p<0.01'); ylim([-3,3.5])
for i=1:length(topbot_tvals) ; if cardob_pvals(i) < 0.01; text(freqs(i),3.3,'*'); end; end
subplot(2,6,6) ; plot(1,'k') ; hold on ; plot(2,'m'); legend({'cardinal orientations','oblique orientations'}); 

bands = {find(freqs>=40 & freqs<=70),find(freqs>=28 & freqs<=32),find(freqs>=15 & freqs<=25),find(freqs>=8 & freqs<=12)};
bottom = oblique; top = cardinal;
bandcolors = {[0,0,0.5],[0,0.5,0],[0.5,0,0],[0.25,0.15,0.15]};
band_names = {'mid gamma (40-70Hz','low gamma (28-32Hz)','beta (15-25Hz)','alpha (8-12Hz)'}; 
subplots = [21,22,23,24]; 
for i=1:length(bands)
    subplot(3,10,subplots(i)); 
    b1 = bar(1,squeeze(mean(mean(mean(mean(compdb(:,[2,4],bands{i},bottom))))))); hold on; 
    b1.FaceColor = bandcolors{i}; % b1.LineStyle = ':'; b1.LineWidth = 3; b1.EdgeColor = 'r'; 
    b2 = bar(2,squeeze(mean(mean(mean(mean(compdb(:,[2,4],bands{i},top))))))); 
    b2.FaceColor = bandcolors{i}; 
    
    [h,p,ci,stats] = ttest(squeeze(mean(mean(mean(compdb(:,[2,4],bands{i},bottom),2),3),4)),squeeze(mean(mean(mean(compdb(:,[2,4],bands{i},top),2),3),4))); 
    
    mean_ebar = [squeeze(mean(mean(mean(mean(compdb(:,[2,4],bands{i},bottom)))))),squeeze(mean(mean(mean(mean(compdb(:,[2,4],bands{i},top))))))]; 
    std_ebar = [squeeze(std(mean(mean(mean(compdb(:,[2,4],bands{i},bottom),2),3),4),0,1))/sqrt(11),squeeze(std(mean(mean(mean(compdb(:,[2,4],bands{i},top),2),3),4),0,1))/sqrt(11)]; 
    errorbar([1,2],mean_ebar,std_ebar,'k.'); if i==1;  ylabel('dB'); end
    set(gca,'XTick',[1,2],'XTickLabel',{'oblique','cardinal'});  xlim([0.5,2.5]); xlabel(band_names{i}); title([format_t(stats.tstat),', ' format_p(p)]);
end
set(gcf,'Position',[100 100 2000 450])


% FIGURE 3
% orientation correlations
retcorrs = corr(squeeze(mean(mean(compdb(:,[1,3],:,:),1),2))'); 
for i=1:size(retcorrs,1); for j=1:size(retcorrs,2) ; if i<j ; retcorrs(i,j) = 0; end ;end; end
subplot(2,6,1) ; imagesc(freqs,freqs,retcorrs,[-1,1]) ; colormap jet; h = colorbar; title(h,'rho'); xlabel('frequency(hz)'); ylabel('frequency(hz)'); 
title('retinotopic tuning corr.');

orientcorrs = corr(squeeze(mean(mean(compdb(:,[2,4],:,:),1),2))'); 
for i=1:size(orientcorrs,1); for j=1:size(orientcorrs,2) ; if i<j ; orientcorrs(i,j) = 0; end ;end; end
subplot(2,6,2) ; imagesc(freqs,freqs,orientcorrs,[-1,1]) ; colormap jet; h = colorbar; title(h,'rho'); xlabel('frequency(hz)'); ylabel('frequency(hz)'); 
title('orientation tuning corr.');

% orientation correlations
clear ocorrs; 
for i=1:11
   ocorrs(i,:,:) = corr(squeeze(mean(compdb(i,[2,4],:,:),2))'); 
   rcorrs(i,:,:) = corr(squeeze(mean(compdb(i,[1,3],:,:),2))'); 
    
end

for i=1:size(rcorrs,2)
    for j=1:size(rcorrs,3)
     [h,p,ci,stats] = ttest(squeeze(rcorrs(:,i,j)),squeeze(ocorrs(:,i,j))); 
     corrts(i,j) = stats.tstat; 
     corrps(i,j) = p;       
    end
end

sigmat = (corrts.*double(corrps<0.05));
for i=1:size(sigmat,1)
    for j=1:size(sigmat,2)
        if j>i ; sigmat(i,j) = NaN; elseif j==i ; sigmat(i,j) = NaN; end
    end
end

subplot(2,6,3);
imagesc(freqs,freqs,sigmat,[-5,5]) ; colormap jet; h = colorbar ; title(h,'t-value'); title('tuning corr. difference'); xlabel('frequency(hz)') ;ylabel('frequency(hz)'); 
set(gcf,'Position',[100 100 2000 450])

%{
subplot(2,6,8); 
bar(1,squeeze(mean(mean(mean(rcorrs(:,freqs>40 & freqs<70,freqs>8 & freqs<=25),2),3)))); hold on ; 
bar(2,squeeze(mean(mean(mean(ocorrs(:,freqs>40 & freqs<70,freqs>8 & freqs<=25),2),3))));
corr_means = [squeeze(mean(mean(mean(rcorrs(:,freqs>40 & freqs<70,freqs>8 & freqs<=25),2),3))),squeeze(mean(mean(mean(ocorrs(:,freqs>40 & freqs<70,freqs>8 & freqs<=25),2),3)))];
corr_stds = [squeeze(std(mean(mean(rcorrs(:,freqs>40 & freqs<70,freqs>8 & freqs<=25),2),3),0,1))/sqrt(11),squeeze(std(mean(mean(ocorrs(:,freqs>40 & freqs<70,freqs>8 & freqs<=25),2),3),0,1))/sqrt(11)]; 
errorbar([1,2],corr_means,corr_stds,'k.'); set(gca,'XTick',[1,2],'XTickLabel',{'retinotopy','orientation'}); xlabel('tuning type'); ylabel('rho'); title('rho alpha/beta vs gamma');
%}
plot(squeeze(mean(mean(mean(abs(compdb(:,[2,4],:,:)),1),2),3)))

wcount = 1;
for w=1:4:301
ideal = zeros(1,360); ideal(vertical) = 1; ideal(horizontal) = -1;
%ideal = imfilter(ideal,fspecial('gaussian',[1,w],w/2));
ideal = smooth(ideal,w)'; 
%figure,subplot(1,2,1) ; plot(ideal); 
otuning = squeeze(mean(mean(compdb(:,[2,4],:,:),1),2)); 
ocorrs = corr(ideal',otuning'); 
sb_otuning = squeeze(mean(compdb(:,[2,4],:,:),2)); 
for i=1:11 ; sb_ocorrs(i,:) = corr(ideal',squeeze(sb_otuning(i,:,:))'); end
%subplot(1,2,2); 
%shadedErrorBar(freqs,squeeze(mean(sb_ocorrs,1)),squeeze(std(sb_ocorrs,0,1))/sqrt(11),{'k'}) ; hline(0,'k'); ylim([-.5,.5]);
gcorrs(wcount) = squeeze(mean(mean(sb_ocorrs(:,1,25:35),1),3)); 
wcount =wcount + 1;
corrfreqs(wcount,:) = squeeze(mean(sb_ocorrs,1)); 
end

%{
subplot(2,3,1); imagesc(1:360,freqs,squeeze(mean(mean(compdb(:,[1,3],:,:),1),2))) ; axis xy ; colormap jet; h=colorbar; title(h,'dB'); xlabel('wedge location (deg)'); ylabel('frequency(hz)'); title('EEG retinotopic tuning');
subplot(2,3,2); imagesc(1:360,freqs,squeeze(mean(mean(compdb(:,[2,4],:,:),1),2))) ; axis xy ; colormap jet; h=colorbar; title(h,'dB'); xlabel('bars orientation (deg)'); ylabel('frequency(hz)'); title('EEG orientation tuning'); 

bandtitles = {'mid gamma (40-70Hz)','low gamma (28-32Hz)','beta (15-25Hz)','alpha (8-12Hz)'};
quadtitles = {'bot right','bot left','top left','top right'};

subplot(2,3,4); 
plot(1:360,squeeze(mean(mean(mean(compdb(:,[1,3],20:35,:),1),2),3)),'Color',[0,0,0.5],'LineWidth',1); hold on; 
plot(1:360,squeeze(mean(mean(mean(compdb(:,[1,3],14:15,:),1),2),3)),'Color',[0,0.5,0],'LineWidth',1); 
plot(1:360,squeeze(mean(mean(mean(compdb(:,[1,3],8:12,:),1),2),3)),'Color',[0.5,0,0],'LineWidth',1); 
plot(1:360,squeeze(mean(mean(mean(compdb(:,[1,3],4:6,:),1),2),3)),'Color',[0.25,0.15,0.15],'LineWidth',1); hline(0,'k'); xlim([0,360]); 
title('frequency specific retinotopic tuning'); xlabel('wedge location (deg)'); ylabel('dB'); 
for i=1:length(quads) ; vline(quads(i),'r'); text(quads(i),0.8,quadtitles{i}); end

subplot(2,3,5); 
plot(1:360,squeeze(mean(mean(mean(compdb(:,[2,4],20:35,:),1),2),3)),'Color',[0,0,0.5],'LineWidth',1); hold on; 
plot(1:360,squeeze(mean(mean(mean(compdb(:,[2,4],14:15,:),1),2),3)),'Color',[0,0.5,0],'LineWidth',1); 
plot(1:360,squeeze(mean(mean(mean(compdb(:,[2,4],8:12,:),1),2),3)),'Color',[0.5,0,0],'LineWidth',1); 
plot(1:360,squeeze(mean(mean(mean(compdb(:,[2,4],4:6,:),1),2),3)),'Color',[0.25,0.15,0.15],'LineWidth',1); hline(0,'k'); xlim([0,360]); 
legend(bandtitles);title('frequency specific orientation tuning'); xlabel('bars orientation (degree)'); ylabel('dB'); 
for i=1:length(oblique) ; vline(oblique(i),'r'); text(oblique(i),1.9,[num2str(oblique(i)),'\circ']); end


retcorrs = corr(squeeze(mean(mean(compdb(:,[1,3],:,:),1),2))'); 
orientcorrs = corr(squeeze(mean(mean(compdb(:,[2,4],:,:),1),2))'); 
subplot(2,2,1); 
imagesc(retcorrs,[-.6,.6]);
subplot(2,2,2) ; 
imagesc(orientcorrs,[-.6,.6]); 

for i=1:11; subretcorrs(i,:,:) = corr(squeeze(mean(compdb(i,[1,3],:,:),2))'); end
for i=1:11; suborientcorrs(i,:,:) = corr(squeeze(mean(compdb(i,[2,4],:,:),2))'); end

subplot(2,2,3); 
imagesc(squeeze(mean(subretcorrs)),[-.6,.6]);
subplot(2,2,4) ; 
imagesc(squeeze(mean(suborientcorrs)),[-.6,.6]); 
%}