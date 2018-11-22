clear all ; close all ; 
subs = {'lowcont_carmen','lowcont_lyndis','lowcont_mael','lowcont_lorne','lowcont_lyes','lowcont_greg','lowcont_angelina','lowcont_cesar','lowcont_charles','lowcont_alex','lowcont_janie','lowcont_russell','lowcont_anne'} ; 
comps = {[3,4,12,19],[7],[15,21],[4,5,8],[4,11],[5,6],[8,12],[38,11],[5,11,13],[11,20],[8,12],[5,16,21],[2,4,16]} ; 
stims = {'S  2','S  4','S 12','S 14','S 22','S 24','S 32','S 34','S 42','S 44'} ;
bades = {[7,15,26,42],[],[12,19],[],[],[19,56],[],[4,12,19,40],[12,46],[16,19,37],[32],[],[31]} ;

% methods figures: 

% use more frequency bands in the noise function 

%cd e:/saved; 
%allstersp = load('allstersp'); allstersp = allstersp.allstersp; 

%alpha_template = squeeze(mean(mean(mean(allstersp(:,1:6,5:12,:),1),2),3)); 
%gamma_template = squeeze(mean(mean(mean(allstersp(:,1:6,16:35,:),1),2),3)); 

for sub=1:length(subs);
    
    cd(['E:/jly_orientation/',subs{sub}]) ;  

    merged = pop_loadset('merged.set') ;
     
    icaw = load('icaw') ; icaw = icaw.icaw ; winv = pinv(icaw{1}*icaw{2}) ; 
    newmerged = merged ; newmerged.data = icaw{1}*icaw{2}*merged.data ; 
    goodcs = comps{sub} ; 
    badsi = load('badsi'); badsi = badsi.badsi; 
    
    % FIND THE GOOD COMPONENTS
    %{
    clear stersp sorttrials
    for s=1:10 ; disp(s) ; 
        allep = pop_epoch(newmerged,{stims{s}},[-1.5,19.5]) ; 
        
        %figure,plot(squeeze(std(allep.data(:,:,:),0,1)))
        if s==1
        [~,sorttrials(s,:)] = sort(mean(squeeze(std(allep.data(:,:,:),0,1)),1),'descend'); 
        else
        [~,sorttrials(s,:)] = sort(mean(squeeze(std(allep.data(:,:,1:size(sorttrials,2)),0,1)),1),'descend'); 
        end
        
    for i=1:64
        for k=1:size(allep.data,3)
            [stersp(s,i,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,k)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',NaN,'verbose','off','timesout',600) ; 
        end
    end
    end
    bersp = stersp - repmat(mean(stersp(:,:,:,:,times<0),5),[1,1,1,1,600]) ; 
    
    alphabersp = squeeze(mean(bersp(:,:,:,freqs>10 & freqs<22,:),4)); 
    gammabersp = squeeze(mean(bersp(:,:,:,freqs>35 & freqs<70,:),4)); 

    clear alphacorrs gammacorrs; 
    for i=1:size(alphabersp,1)
        for j=1:size(alphabersp,2)
            for k=1:size(alphabersp,3)
                alphacorrs(i,j,k) = corr2(gamma_template,squeeze(gammabersp(i,j,k,:))); 
                gammacorrs(i,j,k) = corr2(alpha_template,squeeze(alphabersp(i,j,k,:))); 
            end
        end
    end
    
    bothcorrs = gammacorrs*2 + -alphacorrs; 
    meanbothcorrs = squeeze(mean(mean(bothcorrs,1),3)); 
    [sv,si] = sort(meanbothcorrs,'descend'); 
    
    %nrf_corr_inds = si; save('nrf_corr_inds','nrf_corr_inds'); 
    
    for i=1:6 ; subplot(2,3,i) ; imagesc(squeeze(mean(mean((stersp(i,:,sorttrials(i,:),20:end,:)),2),4))) ; end

    %}
    % sort the components 
    
    %allbersp(sub,:,:,:) = squeeze(mean(mean(bersp,2),3)); 
 
    
    
    % process the good components
    
    nrf_corr_inds = load('nrf_corr_inds'); nrf_corr_inds = nrf_corr_inds.nrf_corr_inds; 
    clear mstersp 
    for s=1:length(stims) ; disp(s) ; 
        allep = pop_epoch(newmerged,{stims{s}},[-2,19]) ; 
        [sv,si] = sort(mean(std(allep.data,0,2),1),'descend');
        
    for i=1:length(comps{sub})
            [mstersp(s,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(comps{sub}(i),:,si(3:end))),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',60,'baseline',0,'verbose','off','timesout',2000) ; 
    end
    end
    % sort the components 
    figure,for i=1:length(comps{sub}) ;subplot(2,3,i) ; imagesc(squeeze(mean(mstersp(:,i,:,:),1)),[-2,2]);  end; suptitle(subs{sub});

    allstersp(sub,:,:,:) = squeeze(mean(mstersp(:,:,:,:),2)); 
    
    
    
    %{
    time_inds = find(times > 0 & times < 18);    
    gtime_inds = find(times>0.5 & times<18);    

     
    horizontal = [355:359,1:5,175:185];
    vertical = [85:95,265:275];
    oblique = [40:50,130:140,220:230,310:320];
    cardinal = [1:10,85:95,175:185,265:275];
    quads = [315,225,135,45];
    bottom = 225:315; top = 45:135;
    right = [315:359,1:45]; left = [135:225];
    
    s2_angles = [90:-1:0,359:-1:45];
    s4_angles = [45:1:359,0:1:90];

    tbersp = bersp(:,:,:,:,time_inds); 
    clear res_angles
    res_angles(1,:) = imresize(s2_angles,[1,length(tbersp)],'nearest'); 
    res_angles(2,:) = imresize(s4_angles,[1,length(tbersp)],'nearest'); 
    res_angles(3,:) = imresize(s2_angles,[1,length(tbersp)],'nearest'); 
    res_angles(4,:) = imresize(s4_angles,[1,length(tbersp)],'nearest'); 
    res_angles(5,:) = imresize(s2_angles,[1,length(tbersp)],'nearest'); 
    res_angles(6,:) = imresize(s4_angles,[1,length(tbersp)],'nearest'); 
    res_angles(7,:) = imresize(s2_angles,[1,length(tbersp)],'nearest'); 
    res_angles(8,:) = imresize(s4_angles,[1,length(tbersp)],'nearest'); 
    res_angles(9,:) = imresize(s2_angles,[1,length(tbersp)],'nearest'); 
    res_angles(10,:) = imresize(s4_angles,[1,length(tbersp)],'nearest'); 
    
    tbersp = bersp(:,:,:,:,gtime_inds); 
    res_angles = res_angles(:,(length(time_inds)-length(gtime_inds)+1):end); 
    
    tbersp_amp = zeros(size(tbersp,1),size(tbersp,2),size(tbersp,3),size(tbersp,4),359);    
    for i=0:359
        for j=1:10
            bersp_ind = find(i==res_angles(j,:)); 
            tbersp_amp(j,:,:,:,i+1) = squeeze(mean(tbersp(j,:,:,:,bersp_ind),5)); 
        end
    end
    save('tbersp_amp','tbersp_amp'); 
    clear badsi; 
    for i=1:10
        [sv,si] = sort(squeeze(abs(zscore(mean(mean(mean(abs(tbersp_amp(i,:,:,:,:)),4),5),2),[],3))),'descend'); 
        badsi(i,:) = si; 
        sub_db(sub,i,:,:) = squeeze(mean(mean(tbersp_amp(i,:,si(2:end),:,:),3),2)); 
    end
    save('badsi','badsi'); 
    %}
end
figure,imagesc(squeeze(mean(mean(allstersp,1),2)),[-2,2]) ; colormap jet; axis xy 
mersp = squeeze(mean(mean(allstersp,1),2)); 
bands = {[1:3],[5:7],[13:15],[23:30]};
for i=1:length(bands) ; plot(times,squeeze(mean(mersp(bands{i},:),1))) ; vline([0,18]); hold on; hline(0,'k'); end

ts = zeros(1,size(mersp,2)); 
ts(times>0 & times<18) = 1; 

mgamma = squeeze(mean(mersp(20:30,:),1)); 
malpha = squeeze(mean(mersp(5:8,:),1)); 

plot(ts) ; hold on ; plot(malpha); 



%{
plot(freqs(2:12),zt(2:12),'bo','LineWidth',2); hold on; 
plot(freqs(18:end),zt(18:end),'ro','LineWidth',2); 
plot(freqs(2:12),zt(2:12),'k','LineWidth',1); hold on; 
plot(freqs(18:end),zt(18:end),'k','LineWidth',1); 
xlabel('frequency(hz)'); ylabel('ideal neuronal response A.U'); 
title('neuronal response function (NRF)');
%}












%{
oblique = [40:50,130:140];
vertical = 85:95; 
horizontal = [1:5,175:180];
cardinal = [85:95,1:5,175:180]; 

subdb(:,1,:,:) = squeeze(mean(sub_db(:,1:2,:,:),2)); 
subdb(:,2,:,:) = squeeze(mean(sub_db(:,3:4,:,:),2)); 
subdb(:,3,:,:) = squeeze(mean(sub_db(:,5:6,:,:),2)); 
subdb(:,4,:,:) = squeeze(mean(sub_db(:,7:8,:,:),2)); 
subdb(:,5,:,:) = squeeze(mean(sub_db(:,9:10,:,:),2)); 

cd e:/saved; 
allhigh = load('allhigh'); allhigh = allhigh.allhigh; 
alllow = load('alllow'); alllow = alllow.alllow; 
mdb = (subdb(:,:,:,1:180) + subdb(:,:,:,181:end))/2; 

contrasts = {'100%','30%','15%','5%','1%'};
for i=1:5 ; subplot(3,7,i) ; 
    imagesc(1:360,freqs,squeeze(mean(subdb(:,i,:,:),1))) ; axis xy ; colormap jet; if i==1; xlabel('bar orientation (deg\circ)'); ylabel('frequency(hz)');  end; title(contrasts{i});
    h = colorbar ; title(h,'dB'); 
end ;
%subplot(3,7,7) ; imagesc([-3,3]) ;  h = colorbar ; title(h,'dB'); 

sfreqs = {20:35,14:16,7:12,4:6};
freqtitles = {'mid gamma (40-70Hz)','low gamma (28-32Hz)','beta (15-25Hz)','alpha (8-12Hz)'};
freq_colors = {[0,0,0.5],[0,0.5,0],[0.5,0,0],[0.25,0.15,0.15]};
contrast_labs = {'100','30','15','5','1'};

% row 2 - show all amplitude cardinal vs oblique
for i=1:length(freqtitles)
    subplot(3,7,i+7); 
    for j=1:5
    card_amp = squeeze(mean(mean(subdb(:,j,sfreqs{i},cardinal),3),4)); 
    obl_amp = squeeze(mean(mean(subdb(:,j,sfreqs{i},oblique),3),4)); 
    
    b1 = bar(j*3-1,mean(card_amp)); hold on; 
    b1.FaceColor = freq_colors{i}; 
    b2 = bar(j*3,mean(obl_amp)); 
    b2.FaceColor = freq_colors{i}; 
    errorbar(j*3-1,mean(card_amp),std(card_amp,0,1)/sqrt(13),'k.'); 
    errorbar(j*3,mean(obl_amp),std(obl_amp,0,1)/sqrt(13),'k.'); 
    end
    set(gca,'XTick',2:3:14,'XTickLabel',contrast_labs); 
    xlabel('%contrast'); 
end

clear allts
for i=1:5 
    for j=1:length(freqs)
   [h,p,ci,stats] = ttest(squeeze(mean(mean(subdb(:,i,j,cardinal),4),3)),squeeze(mean(mean(subdb(:,i,j,oblique),4),3))); 
   allts(i,j) = stats.tstat; allps(i,j) = p; 
    end
end

% PEAK FREQUENCY ESTIMATION:
% zoom in on gamma, plot all subjects with peak frequency (low and mid)
mspec = squeeze(mean(mean(mean(subdb,1),2),4)); 
contrast_colors = {[0.1,0.1,0],[0.2,0.2,0],[0.3,0.3,0],[0.4,0.4,0],[0.5,0.5,0]}; 
oblique_color = 'r'; vertical_color = 'c'; horizontal_color = 'g'; cardinal_color = 'm'; 

ideal_mid = fspecial('gaussian',[1,15],3); 
ideal_low = fspecial('gaussian',[1,12],1); 

orientation_inds = {1:180,oblique,cardinal,vertical,horizontal};

clear xcorrs lowxcorrs
for i=1:size(mdb,1)
    for j=1:size(mdb,2)
        for k=1:length(orientation_inds)
            xcorrs(i,j,k,:) = xcorr(squeeze(mean(mdb(i,j,15:40,orientation_inds{k}),4)),ideal_mid,25); 
            lowxcorrs(i,j,k,:) = xcorr(squeeze(mean(mdb(i,j,10:20,orientation_inds{k}),4)),ideal_low,12); 
        end
    end
end

clear peak_midgamma peak_lowgamma; 
for i=1:size(mdb,1)
    for j=1:size(mdb,2)
        for k=1:length(orientation_inds)
            peak_midgamma(i,j,k) = find(max(squeeze(xcorrs(i,j,k,:)))==squeeze(xcorrs(i,j,k,:)),1); 
            peak_lowgamma(i,j,k) = find(max(squeeze(lowxcorrs(i,j,k,:)))==squeeze(lowxcorrs(i,j,k,:)),1); 
        end
    end
end
peak_midgamma = (peak_midgamma - 27)*2 + 45;
peak_lowgamma = (peak_lowgamma - 11)*2 + 29; 

subplot(4,3,1); 
for i=1:5
   shadedErrorBar(freqs(12:35),squeeze(mean(mean(subdb(:,i,12:35,:),1),4)),squeeze(std(mean(subdb(:,i,12:35,:),4),0,1))/sqrt(13),{'Color',contrast_colors{i}}); hold on ;
   xlabel('frequency(hz)'); ylabel('dB'); xlim([23,70]);
    
end
title('mean gamma (all orientations)'); 
subplot(4,3,2) ; for i=1:5 ; plot(1,'Color',contrast_colors{i},'LineWidth',5) ; hold on ; end ; legend({'100% contrast','50% contrast','25% contrast','15% contrast','5% contrast'});

subplot(4,6,7); 
for i=1:5
b1 = bar(i,squeeze(mean(peak_midgamma(:,i,1),1))); hold on ; errorbar(i,squeeze(mean(peak_midgamma(:,i,1),1)),squeeze(std(peak_midgamma(:,i,1),0,1))/sqrt(13),'k.'); 
[~,p,~,~] = ttest(squeeze(peak_midgamma(:,1,1)),squeeze(peak_midgamma(:,i,1))); 
barval = squeeze(mean(peak_midgamma(:,i,1),1)) + squeeze(std(peak_midgamma(:,i,1),0,1))/sqrt(13); 
if p<0.01 ; text(i,barval+.1,'**') ; elseif p<0.05 ; text(i,barval+0.1,'*'); end
b1.FaceColor = contrast_colors{i}; 
end
title('mid gamma, *p<0.05 **p<0.01'); xlabel('contrast(%)'); set(gca,'XTick',1:5,'XTickLabel',[100,50,25,15,5]);  ylabel('peak frequency (hz)'); 

subplot(4,6,8); 
for i=1:5
b1 = bar(i,squeeze(mean(peak_lowgamma(:,i,1),1))); hold on ; errorbar(i,squeeze(mean(peak_lowgamma(:,i,1),1)),squeeze(std(peak_lowgamma(:,i,1),0,1))/sqrt(13),'k.'); 
[~,p,~,~] = ttest(squeeze(peak_lowgamma(:,1,1)),squeeze(peak_lowgamma(:,i,1))); 
barval = squeeze(mean(peak_lowgamma(:,i,1),1)) + squeeze(std(peak_lowgamma(:,i,1),0,1))/sqrt(13); 
if p<0.01 ; text(i,barval+.1,'**') ; elseif p<0.05 ; text(i,barval+0.1,'*'); end
b1.FaceColor = contrast_colors{i}; 
end
title('low gamma, *p<0.05'); xlabel('contrast(%)'); set(gca,'XTick',1:5,'XTickLabel',[100,50,25,15,5]);  ylabel('peak frequency (hz)'); 

subplot(4,6,13); 
for i=1:5
b1 = bar(i*3-1,squeeze(mean(peak_midgamma(:,i,2),1))); hold on; errorbar(i*3-1,squeeze(mean(peak_midgamma(:,i,2),1)),squeeze(std(peak_midgamma(:,i,2),0,1))/sqrt(13),'k.'); 
b1.FaceColor = oblique_color;
b2 = bar(i*3,squeeze(mean(peak_midgamma(:,i,3),1))); errorbar(i*3,squeeze(mean(peak_midgamma(:,i,3),1)),squeeze(std(peak_midgamma(:,i,3),0,1))/sqrt(13),'k.'); 
b2.FaceColor = cardinal_color; 
[~,p,~,~] = ttest(squeeze(peak_midgamma(:,i,2)),squeeze(peak_midgamma(:,i,3))); 
barh = squeeze(mean(peak_midgamma(:,i,2),1)) + squeeze(std(peak_midgamma(:,i,2),0,1))/sqrt(13); 
if p<0.1 ; text(i*3-.5,barh+0.1,'*'); elseif p<0.05 ; text(i*3-0.5,barh+0.1,'**'); end
xticks(i) = i*3-.5; 
end
title('*p<0.1, **p<0.05'); 
ylabel('peak frequency (hz)'); set(gca,'XTick',xticks,'XTickLabel',[100,50,25,15,5]); xlabel('contrast(%)'); 

subplot(4,6,14); 
for i=1:5
b1 = bar(i*3-1,squeeze(mean(peak_midgamma(:,i,2),1))); hold on; errorbar(i*3-1,squeeze(mean(peak_midgamma(:,i,2),1)),squeeze(std(peak_midgamma(:,i,2),0,1))/sqrt(13),'k.'); 
b1.FaceColor = oblique_color;
b2 = bar(i*3,squeeze(mean(peak_midgamma(:,i,4),1))); errorbar(i*3,squeeze(mean(peak_midgamma(:,i,4),1)),squeeze(std(peak_midgamma(:,i,4),0,1))/sqrt(13),'k.'); 
b2.FaceColor = vertical_color; 
[~,p,~,~] = ttest(squeeze(peak_midgamma(:,i,2)),squeeze(peak_midgamma(:,i,4))); 
barh = squeeze(mean(peak_midgamma(:,i,2),1)) + squeeze(std(peak_midgamma(:,i,2),0,1))/sqrt(13); 
if p<0.1 ; text(i*3-.5,barh+0.1,'*'); elseif p<0.05 ; text(i*3-0.5,barh+0.1,'**'); end
xticks(i) = i*3-.5; 
end
title('*p<0.1, **p<0.05'); 
ylabel('peak frequency (hz)'); set(gca,'XTick',xticks,'XTickLabel',[100,50,25,15,5]); xlabel('contrast(%)'); 


subplot(4,6,15); 
for i=1:5
b1 = bar(i*3-1,squeeze(mean(peak_midgamma(:,i,4),1))); hold on; errorbar(i*3-1,squeeze(mean(peak_midgamma(:,i,4),1)),squeeze(std(peak_midgamma(:,i,4),0,1))/sqrt(13),'k.'); 
b1.FaceColor = vertical_color;
b2 = bar(i*3,squeeze(mean(peak_midgamma(:,i,5),1))); errorbar(i*3,squeeze(mean(peak_midgamma(:,i,5),1)),squeeze(std(peak_midgamma(:,i,5),0,1))/sqrt(13),'k.'); 
b2.FaceColor = horizontal_color; 
[~,p,~,~] = ttest(squeeze(peak_midgamma(:,i,4)),squeeze(peak_midgamma(:,i,5))); 
barh = squeeze(mean(peak_midgamma(:,i,4),1)) + squeeze(std(peak_midgamma(:,i,4),0,1))/sqrt(13); 
if p<0.1 ; text(i*3-.5,barh+0.1,'*'); elseif p<0.05 ; text(i*3-0.5,barh+0.1,'**'); end
xticks(i) = i*3-.5; 
end
title('*p<0.1, **p<0.05'); 
ylabel('peak frequency (hz)'); set(gca,'XTick',xticks,'XTickLabel',[100,50,25,15,5]); xlabel('contrast(%)'); 

subplot(4,6,19) ; plot(1,oblique_color,'LineWidth',5) ; hold on ; plot(1,cardinal_color,'LineWidth',5); legend({'oblique','cardinal'});
subplot(4,6,20) ; plot(1,oblique_color,'LineWidth',5) ; hold on ; plot(1,vertical_color,'LineWidth',5); legend({'oblique','vertical'});
subplot(4,6,21) ; plot(1,vertical_color,'LineWidth',5) ; hold on ; plot(1,horizontal_color,'LineWidth',5); legend({'vertical','horizontal'});

set(gcf,'Position',[100 100 1700 900])

cd e:/saved; 

plot(squeeze(mean(mean(allhigh(2:end,:,:),1),3))); hold on ; 
plot(squeeze(mean(mean(alllow(2:end,:,:),1),3))) ; vline([40:50,130:140]);


% FMRI RESULTS

boldhigh = squeeze(mean(allhigh(2:end,:,:),1)); 
boldlow = squeeze(mean(alllow(2:end,:,:),1)); 
bold_contrast_labels = {'100','5'};
% figure 1a contrast orientation tuning curves
%figure 1a
subplot(4,3,1);
shadedErrorBar([],squeeze(mean(mean(allhigh(2:end,:,:),1),3)),squeeze(std(mean(allhigh(2:end,:,:),1),0,3))/sqrt(14),{'Color',[.1,.1,0]}); hold on; 
shadedErrorBar([],squeeze(mean(mean(alllow(2:end,:,:),1),3)),squeeze(std(mean(alllow(2:end,:,:),1),0,3))/sqrt(14),{'Color',[0.4,0.4,0]});
xlabel('angle (deg\circ)'); ylabel('BOLD amplitude (a.u)'); grid on; set(gca,'XTick',[45,90,135,180]); title('BOLD orientation tuning *p<0.05');  xlim([0,180]); 
for i=1:180; [~,p(i),~,~] = ttest(boldhigh(i,:),boldlow(i,:)); if p(i)<0.05; text(i,3.8,'*'); end; end

subplot(4,3,2) ; plot(1,'Color',[.1,.1,0],'LineWidth',5); hold on ; plot(1,'Color',[.4,.4,0],'LineWidth',5); legend({'100% contrast','5% contrast'});


subplot(4,6,7);
b1 = bar(1,mean(mean(boldhigh(oblique,:)))); hold on; errorbar(1,mean(mean(boldhigh(oblique,:))),std(mean(boldhigh(oblique,:),1),0,2)/sqrt(14),'k.'); 
b1.FaceColor = oblique_color; 
b2 = bar(2,mean(mean(boldhigh(cardinal,:)))); errorbar(2,mean(mean(boldhigh(cardinal,:))),std(mean(boldhigh(cardinal,:),1),0,2)/sqrt(14),'k.'); 
b2.FaceColor = cardinal_color;
b3 = bar(4,mean(mean(boldlow(oblique,:)))); errorbar(4,mean(mean(boldlow(oblique,:))),std(mean(boldlow(oblique,:),1),0,2)/sqrt(14),'k.'); 
b3.FaceColor = oblique_color;
b4 = bar(5,mean(mean(boldlow(cardinal,:)))); errorbar(5,mean(mean(boldlow(cardinal,:))),std(mean(boldlow(cardinal,:),1),0,2)/sqrt(14),'k.'); 
b4.FaceColor = cardinal_color;
ymin = double(min([b1.YData, b2.YData, b3.YData, b4.YData])); ymax = double(max([b1.YData, b2.YData, b3.YData, b4.YData])); 
[~,p,~,~] = ttest(mean(boldhigh(oblique,:)),mean(boldhigh(cardinal,:))); if p<0.05 ; formp = ['*',format_p(p)]; else formp = format_p(p); end ;text(1.5,ymax+abs(ymax/10),formp);
[~,p,~,~] = ttest(mean(boldlow(oblique,:)),mean(boldlow(cardinal,:))); if p<0.05 ; formp = ['*',format_p(p)]; else formp = format_p(p); end ; text(4.5,ymax+abs(ymax/10),formp);
set(gca,'XTick',[1,4],'XTickLabel',bold_contrast_labels); xlabel('contrast(%)'); ylabel('BOLD amplitude (AU)'); title('oblique vs cardinal');

subplot(4,6,8); 
b1 = bar(1,mean(mean(boldhigh(oblique,:)))); hold on; errorbar(1,mean(mean(boldhigh(oblique,:))),std(mean(boldhigh(oblique,:),1),0,2)/sqrt(14),'k.'); 
b1.FaceColor = oblique_color; 
b2 = bar(2,mean(mean(boldhigh(vertical,:)))); errorbar(2,mean(mean(boldhigh(vertical,:))),std(mean(boldhigh(vertical,:),1),0,2)/sqrt(14),'k.'); 
b2.FaceColor = vertical_color; 
b3 = bar(4,mean(mean(boldlow(oblique,:)))); errorbar(4,mean(mean(boldlow(oblique,:))),std(mean(boldlow(oblique,:),1),0,2)/sqrt(14),'k.'); 
b3.FaceColor = oblique_color;
b4 = bar(5,mean(mean(boldlow(vertical,:)))); errorbar(5,mean(mean(boldlow(vertical,:))),std(mean(boldlow(vertical,:),1),0,2)/sqrt(14),'k.'); 
b4.FaceColor = vertical_color; 
ymin = double(min([b1.YData, b2.YData, b3.YData, b4.YData])); ymax = double(max([b1.YData, b2.YData, b3.YData, b4.YData])); 
[~,p,~,~] = ttest(mean(boldhigh(oblique,:)),mean(boldhigh(vertical,:))); if p<0.05 ; formp = ['*',format_p(p)]; else formp = format_p(p); end ;text(1.5,ymax+abs(ymax/10),formp);
[~,p,~,~] = ttest(mean(boldlow(oblique,:)),mean(boldlow(vertical,:))); if p<0.05 ; formp = ['*',format_p(p)]; else formp = format_p(p); end ; text(4.5,ymax+abs(ymax/10),formp);
set(gca,'XTick',[1,4],'XTickLabel',bold_contrast_labels); xlabel('contrast(%)'); ylabel('BOLD amplitude (AU)'); title('oblique vs vertical');

subplot(4,6,9); 
b1 = bar(1,mean(mean(boldhigh(vertical,:)))); hold on; errorbar(1,mean(mean(boldhigh(vertical,:))),std(mean(boldhigh(vertical,:),1),0,2)/sqrt(14),'k.'); 
b1.FaceColor = vertical_color; 
b2 = bar(2,mean(mean(boldhigh(horizontal,:)))); errorbar(2,mean(mean(boldhigh(horizontal,:))),std(mean(boldhigh(horizontal,:),1),0,2)/sqrt(14),'k.'); 
b2.FaceColor = horizontal_color; 
b3 = bar(4,mean(mean(boldlow(vertical,:)))); errorbar(4,mean(mean(boldlow(vertical,:))),std(mean(boldlow(vertical,:),1),0,2)/sqrt(14),'k.'); 
b3.FaceColor = vertical_color;
b4 = bar(5,mean(mean(boldlow(horizontal,:)))); errorbar(5,mean(mean(boldlow(horizontal,:))),std(mean(boldlow(horizontal,:),1),0,2)/sqrt(14),'k.'); 
b4.FaceColor = horizontal_color; 
ymin = double(min([b1.YData, b2.YData, b3.YData, b4.YData])); ymax = double(max([b1.YData, b2.YData, b3.YData, b4.YData])); 
[~,p,~,~] = ttest(mean(boldhigh(vertical,:)),mean(boldhigh(horizontal,:))); if p<0.05 ; formp = ['*',format_p(p)]; else formp = format_p(p); end ;text(1.5,ymax+abs(ymax/10),formp);
[~,p,~,~] = ttest(mean(boldlow(vertical,:)),mean(boldlow(horizontal,:))); if p<0.05 ; formp = ['*',format_p(p)]; else formp = format_p(p); end ; text(4.5,ymax+abs(ymax/10),formp);
set(gca,'XTick',[1,4],'XTickLabel',bold_contrast_labels); xlabel('contrast(%)'); ylabel('BOLD amplitude (AU)'); title('vertical vs horizontal');

subplot(4,6,13) ; plot(1,oblique_color,'LineWidth',5); hold on ; plot(1,cardinal_color,'LineWidth',5); legend({'oblique','cardinal'});
subplot(4,6,14) ; plot(1,oblique_color,'LineWidth',5); hold on ; plot(1,vertical_color,'LineWidth',5); legend({'oblique','vertical'});
subplot(4,6,15) ; plot(1,vertical_color,'LineWidth',5); hold on ; plot(1,horizontal_color,'LineWidth',5); legend({'vertical','horizontal'});
set(gcf,'Position',[100 100 1700 900])

% END FIGURE 1



%%% figure 2: eeg orientation tuning
eeg_contrast_names = {'100','50','25','15','5'};
for i=1:5 
    subplot(4,6,i) ;
    imagesc(squeeze(mean(mdb(:,i,:,:),1)),[-3.5,2.5]); axis xy ; colormap jet; 
    set(gca,'XTick',[45,90,135]);  title([eeg_contrast_names{i},'% contrast']);
    if i==1 ; xlabel('angle(deg\circ)'); ylabel('frequency(hz)'); end
end 
subplot(4,6,6) ; imagesc([-3.5,2.5]); h = colorbar ; title(h,'dB'); 

gamma = [find(freqs>= 28 & freqs<=32), find(freqs>=40 & freqs<=70)];
alphabeta = find(freqs>=10& freqs <=15); 

subplot(4,3,4); 
shadedErrorBar(1:180,squeeze(mean(mean(mdb(:,1,alphabeta,:),1),3)),squeeze(std(mean(mdb(:,1,alphabeta,:),3),0,1))/sqrt(13),{'Color',contrast_colors{1}}); hold on ; 
shadedErrorBar(1:180,squeeze(mean(mean(mdb(:,5,alphabeta,:),1),3)),squeeze(std(mean(mdb(:,5,alphabeta,:),3),0,1))/sqrt(13),{'Color',contrast_colors{5}}); hold on ; 
set(gca,'XTick',[45,90,135]) ; xlabel('orientation (deg\circ)'); ylabel('dB'); grid on; title('alpha/beta orientation tuning, *p<0.05'); ylim([-3.5,-1.4]);
for i=1:180 ; [~,p,~,~] = ttest(squeeze(mean(mdb(:,1,alphabeta,i),3)),squeeze(mean(mdb(:,5,alphabeta,i),3))); if p<0.05 ; text(i,-1.5,'*') ; end; end ; xlim([0,180]);

subplot(4,3,5); 
shadedErrorBar(1:180,squeeze(mean(mean(mdb(:,1,gamma,:),1),3)),squeeze(std(mean(mdb(:,1,gamma,:),3),0,1))/sqrt(13),{'Color',contrast_colors{1}}); hold on ; 
shadedErrorBar(1:180,squeeze(mean(mean(mdb(:,5,gamma,:),1),3)),squeeze(std(mean(mdb(:,5,gamma,:),3),0,1))/sqrt(13),{'Color',contrast_colors{5}}); hold on ; 
set(gca,'XTick',[45,90,135]) ; xlabel('orientation (deg\circ)'); ylabel('dB'); grid on; title('gamma orientation tuning, *p<0.05'); ylim([-.5,1.9]);
for i=1:180 ; [~,p,~,~] = ttest(squeeze(mean(mdb(:,1,gamma,i),3)),squeeze(mean(mdb(:,5,gamma,i),3))); if p<0.05 ; text(i,1.8,'*') ; end; end;xlim([0,180]);
subplot(4,3,6) ; plot(1,'Color',contrast_colors{1},'LineWidth',5) ; hold on ; plot(1,'Color',contrast_colors{5},'LineWidth',5) ; legend({'100% contrast','5% contrast'});

subplot(4,6,13); 
b1 = bar(1,squeeze(mean(mean(mean(mdb(:,1,gamma,oblique),1),3),4))); hold on ; errorbar(1,squeeze(mean(mean(mean(mdb(:,1,gamma,oblique),1),3),4)),squeeze(std(mean(mean(mdb(:,1,gamma,oblique),3),4),0,1))/sqrt(13),'k.');
b1.FaceColor = oblique_color; 
b2 = bar(2,squeeze(mean(mean(mean(mdb(:,1,gamma,cardinal),1),3),4))); errorbar(2,squeeze(mean(mean(mean(mdb(:,1,gamma,cardinal),1),3),4)),squeeze(std(mean(mean(mdb(:,1,gamma,cardinal),3),4),0,1))/sqrt(13),'k.');
b2.FaceColor = cardinal_color; 
[~,p,~,~] = ttest(squeeze(mean(mean(mdb(:,1,gamma,oblique),3),4)),squeeze(mean(mean(mdb(:,1,gamma,cardinal),3),4))); if p < 0.05 ; formp = ['*',format_p(p)] ; else formp = format_p(p); end 
ymax = double(max([b1.YData,b2.YData])); text(1,ymax+0.1,formp); 

b1 = bar(4,squeeze(mean(mean(mean(mdb(:,5,gamma,oblique),1),3),4))); hold on ; errorbar(4,squeeze(mean(mean(mean(mdb(:,5,gamma,oblique),1),3),4)),squeeze(std(mean(mean(mdb(:,5,gamma,oblique),3),4),0,1))/sqrt(13),'k.');
b1.FaceColor = oblique_color; 
b2 = bar(5,squeeze(mean(mean(mean(mdb(:,5,gamma,cardinal),1),3),4))); errorbar(5,squeeze(mean(mean(mean(mdb(:,5,gamma,cardinal),1),3),4)),squeeze(std(mean(mean(mdb(:,5,gamma,cardinal),3),4),0,1))/sqrt(13),'k.');
b2.FaceColor = cardinal_color; 
[~,p,~,~] = ttest(squeeze(mean(mean(mdb(:,5,gamma,oblique),3),4)),squeeze(mean(mean(mdb(:,5,gamma,cardinal),3),4))); if p < 0.05 ; formp = ['*',format_p(p)] ; else formp = format_p(p); end 
ymax = double(max([b1.YData,b2.YData])); text(4,ymax+0.25,formp); 
title('oblique vs cardinal (gamma)'); set(gca,'XTick',[1.5,4.5],'XTickLabel',{'100','5'}) ; xlabel('contrast(%)'); ylabel('dB'); 

subplot(4,6,14); 
b1 = bar(1,squeeze(mean(mean(mean(mdb(:,1,gamma,oblique),1),3),4))); hold on ; errorbar(1,squeeze(mean(mean(mean(mdb(:,1,gamma,oblique),1),3),4)),squeeze(std(mean(mean(mdb(:,1,gamma,oblique),3),4),0,1))/sqrt(13),'k.');
b1.FaceColor = oblique_color; 
b2 = bar(2,squeeze(mean(mean(mean(mdb(:,1,gamma,vertical),1),3),4))); errorbar(2,squeeze(mean(mean(mean(mdb(:,1,gamma,vertical),1),3),4)),squeeze(std(mean(mean(mdb(:,1,gamma,vertical),3),4),0,1))/sqrt(13),'k.');
b2.FaceColor = vertical_color; 
[~,p,~,~] = ttest(squeeze(mean(mean(mdb(:,1,gamma,oblique),3),4)),squeeze(mean(mean(mdb(:,1,gamma,vertical),3),4))); if p < 0.05 ; formp = ['*',format_p(p)] ; else formp = format_p(p); end 
ymax = double(max([b1.YData,b2.YData])); text(1,ymax+0.1,formp); 

b1 = bar(4,squeeze(mean(mean(mean(mdb(:,5,gamma,oblique),1),3),4))); hold on ; errorbar(4,squeeze(mean(mean(mean(mdb(:,5,gamma,oblique),1),3),4)),squeeze(std(mean(mean(mdb(:,5,gamma,oblique),3),4),0,1))/sqrt(13),'k.');
b1.FaceColor = oblique_color; 
b2 = bar(5,squeeze(mean(mean(mean(mdb(:,5,gamma,vertical),1),3),4))); errorbar(5,squeeze(mean(mean(mean(mdb(:,5,gamma,vertical),1),3),4)),squeeze(std(mean(mean(mdb(:,5,gamma,vertical),3),4),0,1))/sqrt(13),'k.');
b2.FaceColor = vertical_color; 
[~,p,~,~] = ttest(squeeze(mean(mean(mdb(:,5,gamma,oblique),3),4)),squeeze(mean(mean(mdb(:,5,gamma,vertical),3),4))); if p < 0.05 ; formp = ['*',format_p(p)] ; else formp = format_p(p); end 
ymax = double(max([b1.YData,b2.YData])); text(4,ymax+0.25,formp); 
title('oblique vs vertical (gamma)'); set(gca,'XTick',[1.5,4.5],'XTickLabel',{'100','5'}) ; xlabel('contrast(%)'); ylabel('dB'); 

subplot(4,6,15); 
b1 = bar(1,squeeze(mean(mean(mean(mdb(:,1,gamma,vertical),1),3),4))); hold on ; errorbar(1,squeeze(mean(mean(mean(mdb(:,1,gamma,vertical),1),3),4)),squeeze(std(mean(mean(mdb(:,1,gamma,vertical),3),4),0,1))/sqrt(13),'k.');
b1.FaceColor = vertical_color; 
b2 = bar(2,squeeze(mean(mean(mean(mdb(:,1,gamma,horizontal),1),3),4))); errorbar(2,squeeze(mean(mean(mean(mdb(:,1,gamma,horizontal),1),3),4)),squeeze(std(mean(mean(mdb(:,1,gamma,horizontal),3),4),0,1))/sqrt(13),'k.');
b2.FaceColor = horizontal_color; 
[~,p,~,~] = ttest(squeeze(mean(mean(mdb(:,1,gamma,vertical),3),4)),squeeze(mean(mean(mdb(:,1,gamma,horizontal),3),4))); if p < 0.05 ; formp = ['*',format_p(p)] ; else formp = format_p(p); end 
ymax = double(max([b1.YData,b2.YData])); text(1,ymax+0.1,formp); 

b1 = bar(4,squeeze(mean(mean(mean(mdb(:,5,gamma,vertical),1),3),4))); hold on ; errorbar(4,squeeze(mean(mean(mean(mdb(:,5,gamma,vertical),1),3),4)),squeeze(std(mean(mean(mdb(:,5,gamma,vertical),3),4),0,1))/sqrt(13),'k.');
b1.FaceColor = vertical_color; 
b2 = bar(5,squeeze(mean(mean(mean(mdb(:,5,gamma,horizontal),1),3),4))); errorbar(5,squeeze(mean(mean(mean(mdb(:,5,gamma,horizontal),1),3),4)),squeeze(std(mean(mean(mdb(:,5,gamma,horizontal),3),4),0,1))/sqrt(13),'k.');
b2.FaceColor = horizontal_color; 
[~,p,~,~] = ttest(squeeze(mean(mean(mdb(:,5,gamma,vertical),3),4)),squeeze(mean(mean(mdb(:,5,gamma,horizontal),3),4))); if p < 0.05 ; formp = ['*',format_p(p)] ; else formp = format_p(p); end 
ymax = double(max([b1.YData,b2.YData])); text(4,ymax+0.25,formp); 
title('vertical vs horizontal (gamma)'); set(gca,'XTick',[1.5,4.5],'XTickLabel',{'100','5'}) ; xlabel('contrast(%)'); ylabel('dB'); 



subplot(4,6,19); 
b1 = bar(1,squeeze(mean(mean(mean(mdb(:,1,alphabeta,oblique),1),3),4))); hold on ; errorbar(1,squeeze(mean(mean(mean(mdb(:,1,alphabeta,oblique),1),3),4)),squeeze(std(mean(mean(mdb(:,1,alphabeta,oblique),3),4),0,1))/sqrt(13),'k.');
b1.FaceColor = oblique_color; 
b2 = bar(2,squeeze(mean(mean(mean(mdb(:,1,alphabeta,cardinal),1),3),4))); errorbar(2,squeeze(mean(mean(mean(mdb(:,1,alphabeta,cardinal),1),3),4)),squeeze(std(mean(mean(mdb(:,1,alphabeta,cardinal),3),4),0,1))/sqrt(13),'k.');
b2.FaceColor = cardinal_color; 
[~,p,~,~] = ttest(squeeze(mean(mean(mdb(:,1,alphabeta,oblique),3),4)),squeeze(mean(mean(mdb(:,1,alphabeta,cardinal),3),4))); if p < 0.05 ; formp = ['*',format_p(p)] ; else formp = format_p(p); end 
ymin = double(min([b1.YData,b2.YData])); text(1,ymin-0.3,formp); 

b1 = bar(4,squeeze(mean(mean(mean(mdb(:,5,alphabeta,oblique),1),3),4))); hold on ; errorbar(4,squeeze(mean(mean(mean(mdb(:,5,alphabeta,oblique),1),3),4)),squeeze(std(mean(mean(mdb(:,5,alphabeta,oblique),3),4),0,1))/sqrt(13),'k.');
b1.FaceColor = oblique_color; 
b2 = bar(5,squeeze(mean(mean(mean(mdb(:,5,alphabeta,cardinal),1),3),4))); errorbar(5,squeeze(mean(mean(mean(mdb(:,5,alphabeta,cardinal),1),3),4)),squeeze(std(mean(mean(mdb(:,5,alphabeta,cardinal),3),4),0,1))/sqrt(13),'k.');
b2.FaceColor = cardinal_color; 
[~,p,~,~] = ttest(squeeze(mean(mean(mdb(:,5,alphabeta,oblique),3),4)),squeeze(mean(mean(mdb(:,5,alphabeta,cardinal),3),4))); if p < 0.05 ; formp = ['*',format_p(p)] ; else formp = format_p(p); end 
ymin = double(min([b1.YData,b2.YData])); text(4,ymin-0.3,formp); 
title('oblique vs cardinal (alpha/beta)'); set(gca,'XTick',[1.5,4.5],'XTickLabel',{'100','5'}) ; xlabel('contrast(%)'); ylabel('dB'); 


subplot(4,6,20); 
b1 = bar(1,squeeze(mean(mean(mean(mdb(:,1,alphabeta,oblique),1),3),4))); hold on ; errorbar(1,squeeze(mean(mean(mean(mdb(:,1,alphabeta,oblique),1),3),4)),squeeze(std(mean(mean(mdb(:,1,alphabeta,oblique),3),4),0,1))/sqrt(13),'k.');
b1.FaceColor = oblique_color; 
b2 = bar(2,squeeze(mean(mean(mean(mdb(:,1,alphabeta,vertical),1),3),4))); errorbar(2,squeeze(mean(mean(mean(mdb(:,1,alphabeta,vertical),1),3),4)),squeeze(std(mean(mean(mdb(:,1,alphabeta,vertical),3),4),0,1))/sqrt(13),'k.');
b2.FaceColor = vertical_color; 
[~,p,~,~] = ttest(squeeze(mean(mean(mdb(:,1,alphabeta,oblique),3),4)),squeeze(mean(mean(mdb(:,1,alphabeta,vertical),3),4))); if p < 0.05 ; formp = ['*',format_p(p)] ; else formp = format_p(p); end 
ymin = double(min([b1.YData,b2.YData])); text(1,ymin-0.3,formp); 

b1 = bar(4,squeeze(mean(mean(mean(mdb(:,5,alphabeta,oblique),1),3),4))); hold on ; errorbar(4,squeeze(mean(mean(mean(mdb(:,5,alphabeta,oblique),1),3),4)),squeeze(std(mean(mean(mdb(:,5,alphabeta,oblique),3),4),0,1))/sqrt(13),'k.');
b1.FaceColor = oblique_color; 
b2 = bar(5,squeeze(mean(mean(mean(mdb(:,5,alphabeta,vertical),1),3),4))); errorbar(5,squeeze(mean(mean(mean(mdb(:,5,alphabeta,vertical),1),3),4)),squeeze(std(mean(mean(mdb(:,5,alphabeta,vertical),3),4),0,1))/sqrt(13),'k.');
b2.FaceColor = vertical_color; 
[~,p,~,~] = ttest(squeeze(mean(mean(mdb(:,5,alphabeta,oblique),3),4)),squeeze(mean(mean(mdb(:,5,alphabeta,vertical),3),4))); if p < 0.05 ; formp = ['*',format_p(p)] ; else formp = format_p(p); end 
ymin = double(min([b1.YData,b2.YData])); text(4,ymin-0.3,formp); 
title('oblique vs vertical (alpha/beta)'); set(gca,'XTick',[1.5,4.5],'XTickLabel',{'100','5'}) ; xlabel('contrast(%)'); ylabel('dB'); 

subplot(4,6,21); 
b1 = bar(1,squeeze(mean(mean(mean(mdb(:,1,alphabeta,vertical),1),3),4))); hold on ; errorbar(1,squeeze(mean(mean(mean(mdb(:,1,alphabeta,vertical),1),3),4)),squeeze(std(mean(mean(mdb(:,1,alphabeta,vertical),3),4),0,1))/sqrt(13),'k.');
b1.FaceColor = vertical_color; 
b2 = bar(2,squeeze(mean(mean(mean(mdb(:,1,alphabeta,horizontal),1),3),4))); errorbar(2,squeeze(mean(mean(mean(mdb(:,1,alphabeta,horizontal),1),3),4)),squeeze(std(mean(mean(mdb(:,1,alphabeta,horizontal),3),4),0,1))/sqrt(13),'k.');
b2.FaceColor = horizontal_color; 
[~,p,~,~] = ttest(squeeze(mean(mean(mdb(:,1,alphabeta,vertical),3),4)),squeeze(mean(mean(mdb(:,1,alphabeta,horizontal),3),4))); if p < 0.05 ; formp = ['*',format_p(p)] ; else formp = format_p(p); end 
ymax = double(max([b1.YData,b2.YData])); text(1,ymax-0.25,formp); 

b1 = bar(4,squeeze(mean(mean(mean(mdb(:,5,alphabeta,vertical),1),3),4))); hold on ; errorbar(4,squeeze(mean(mean(mean(mdb(:,5,alphabeta,vertical),1),3),4)),squeeze(std(mean(mean(mdb(:,5,alphabeta,vertical),3),4),0,1))/sqrt(13),'k.');
b1.FaceColor = vertical_color; 
b2 = bar(5,squeeze(mean(mean(mean(mdb(:,5,alphabeta,horizontal),1),3),4))); errorbar(5,squeeze(mean(mean(mean(mdb(:,5,alphabeta,horizontal),1),3),4)),squeeze(std(mean(mean(mdb(:,5,alphabeta,horizontal),3),4),0,1))/sqrt(13),'k.');
b2.FaceColor = horizontal_color; 
[~,p,~,~] = ttest(squeeze(mean(mean(mdb(:,5,alphabeta,vertical),3),4)),squeeze(mean(mean(mdb(:,5,alphabeta,horizontal),3),4))); if p < 0.05 ; formp = ['*',format_p(p)] ; else formp = format_p(p); end 
ymin = double(min([b1.YData,b2.YData])); text(4,ymin-0.25,formp); 
title('vertical vs horizontal (alphabeta)'); set(gca,'XTick',[1.5,4.5],'XTickLabel',{'100','5'}) ; xlabel('contrast(%)'); ylabel('dB'); 

subplot(4,6,22) ; plot(1,oblique_color,'LineWidth',5); hold on ; plot(1,cardinal_color,'LineWidth',5); legend({'oblique','cardinal'});
subplot(4,6,23) ; plot(1,oblique_color,'LineWidth',5); hold on ; plot(1,vertical_color,'LineWidth',5); legend({'oblique','vertical'});
subplot(4,6,24) ; plot(1,vertical_color,'LineWidth',5); hold on ; plot(1,horizontal_color,'LineWidth',5); legend({'vertical','horizontal'});


set(gcf,'Position',[100 100 1700 900])


% figure 4 (comparison) 
figure
clear highcorrs lowcorrs highps lowps
for i=1:5 
    for j=1:60
        [highcorrs(i,j),highps(i,j)] = corr(mean(boldhigh,2),squeeze(mean(mdb(:,i,j,:),1))); 
        [lowcorrs(i,j),lowps(i,j)] = corr(mean(boldlow,2),squeeze(mean(mdb(:,i,j,:),1))); 
    end
end

subplot(4,3,1) ; 
imagesc(freqs,1:5,highcorrs,[-1,1]) ; colormap jet; xlabel('frequency(hz)') ; ylabel('EEG contrast(%)') ;set(gca,'YTick',1:5,'YTickLabel',[100,50,25,15,5]); hline([0.5:5.5],'k'); 
title('EEG (all contrast) vs BOLD (100% contrast)'); 
subplot(4,3,2); 
imagesc(freqs,1:5,lowcorrs,[-1,1]) ; colormap jet; xlabel('frequency(hz)') ; ylabel('EEG contrast(%)') ;set(gca,'YTick',1:5,'YTickLabel',[100,50,25,15,5]); hline([0.5:5.5],'k'); 
title('EEG (all contrast) vs BOLD (5% contrast)'); 
subplot(4,3,3) ; imagesc([-1,1]); h = colorbar; title(h,'rho'); 

subplot(4,3,4) ; 
plot(freqs,squeeze(highcorrs(1,:)),'Color',contrast_colors{1},'LineWidth',2); hold on ; hline(0,'k'); 
plot(freqs,squeeze(lowcorrs(5,:)),'Color',contrast_colors{5},'LineWidth',2); hold on ; hline(0,'k'); 
xlim([0,120]); ylabel('rho'); xlabel('frequency(hz)'); title('EEG-BOLD tuning correlation'); 

subplot(4,3,5) ; 
plot(freqs,squeeze(highcorrs(1,:)).^2,'Color',contrast_colors{1},'LineWidth',2); hold on ; hline(0,'k'); 
plot(freqs,squeeze(lowcorrs(5,:)).^2,'Color',contrast_colors{5},'LineWidth',2); hold on ; hline(0,'k'); 
xlim([0,120]); ylabel('rho^{2}'); xlabel('frequency(hz)'); title('BOLD tuning variance explained'); 

subplot(4,3,6); plot(1,'Color',contrast_colors{1},'LineWidth',5); hold on; plot(1,'Color',contrast_colors{5},'LineWidth',5); legend({'100% contrast','5% contrast'});

orientations = zeros(1,180) ; orientations([oblique,vertical,horizontal]) = 1; other_orientations = find(orientations==0); 

subplot(4,6,13); 
plot(squeeze(mean(mean(mdb(:,1,gamma,other_orientations),1),3)),squeeze(mean(boldhigh(other_orientations,:),2)),'d','LineWidth',1,'Color',[0.3,0.3,0.3]);  hold on ; 
plot(squeeze(mean(mean(mdb(:,1,gamma,oblique),1),3)),squeeze(mean(boldhigh(oblique,:),2)),'d','LineWidth',2,'Color',oblique_color);
plot(squeeze(mean(mean(mdb(:,1,gamma,vertical),1),3)),squeeze(mean(boldhigh(vertical,:),2)),'d','LineWidth',2,'Color',vertical_color); 
plot(squeeze(mean(mean(mdb(:,1,gamma,horizontal),1),3)),squeeze(mean(boldhigh(horizontal,:),2)),'d','LineWidth',2,'Color',horizontal_color); 
x = double(squeeze(mean(mean(mdb(:,1,gamma,:),1),3))); y = double(squeeze(mean(boldhigh,2))); 
B = [ones(size(x(:)))  x(:)]\y(:);
y_fit = [ones(size(x(:)))  x(:)]*B;
plot(x,y_fit,'k','LineWidth',3); xlabel('gamma 100% contrast (dB)'); ylabel('BOLD 100% contrast (AU)'); title('100%contrast (gamma vs BOLD)'); 
[rho,p] = corr(x,y); text(min(x)-abs(min(x))/10,max(y)+abs(max(y))/10,format_rho(rho)); xlim([min(x)-abs(min(x))/10,max(x)+abs(max(x))/10]); ylim([min(y)-abs(min(y))/10,max(y)+abs(max(y))/4]); 

subplot(4,6,14); 
plot(squeeze(mean(mean(mdb(:,1,alphabeta,other_orientations),1),3)),squeeze(mean(boldhigh(other_orientations,:),2)),'d','LineWidth',1,'Color',[0.3,0.3,0.3]);  hold on ; 
plot(squeeze(mean(mean(mdb(:,1,alphabeta,oblique),1),3)),squeeze(mean(boldhigh(oblique,:),2)),'d','LineWidth',2,'Color',oblique_color);
plot(squeeze(mean(mean(mdb(:,1,alphabeta,vertical),1),3)),squeeze(mean(boldhigh(vertical,:),2)),'d','LineWidth',2,'Color',vertical_color); 
plot(squeeze(mean(mean(mdb(:,1,alphabeta,horizontal),1),3)),squeeze(mean(boldhigh(horizontal,:),2)),'d','LineWidth',2,'Color',horizontal_color); 
x = double(squeeze(mean(mean(mdb(:,1,alphabeta,:),1),3))); y = double(squeeze(mean(boldhigh,2))); 
B = [ones(size(x(:)))  x(:)]\y(:);
y_fit = [ones(size(x(:)))  x(:)]*B;
plot(x,y_fit,'k','LineWidth',3); xlabel('alpha/beta 100% contrast (dB)'); ylabel('BOLD 100% contrast (AU)'); title('100%contrast (alpha/beta vs BOLD)'); 
[rho,p] = corr(x,y); text(min(x)-abs(min(x))/10,max(y)+abs(max(y))/10,format_rho(rho));   xlim([min(x)-abs(min(x))/10,max(x)+abs(max(x))/10]); ylim([min(y)-abs(min(y))/10,max(y)+abs(max(y))/4]); 

subplot(4,6,15); 
plot(squeeze(mean(mean(mdb(:,5,gamma,other_orientations),1),3)),squeeze(mean(boldlow(other_orientations,:),2)),'d','LineWidth',1,'Color',[0.3,0.3,0.3]);  hold on ; 
plot(squeeze(mean(mean(mdb(:,5,gamma,oblique),1),3)),squeeze(mean(boldlow(oblique,:),2)),'d','LineWidth',2,'Color',oblique_color);
plot(squeeze(mean(mean(mdb(:,5,gamma,vertical),1),3)),squeeze(mean(boldlow(vertical,:),2)),'d','LineWidth',2,'Color',vertical_color); 
plot(squeeze(mean(mean(mdb(:,5,gamma,horizontal),1),3)),squeeze(mean(boldlow(horizontal,:),2)),'d','LineWidth',2,'Color',horizontal_color); 
x = double(squeeze(mean(mean(mdb(:,5,gamma,:),1),3))); y = double(squeeze(mean(boldlow,2))); 
B = [ones(size(x(:)))  x(:)]\y(:);
y_fit = [ones(size(x(:)))  x(:)]*B;
plot(x,y_fit,'k','LineWidth',3); xlabel('gamma 5% contrast (dB)'); ylabel('BOLD 5% contrast (AU)'); title('5%contrast (gamma vs BOLD)'); 
[rho,p] = corr(x,y); text(min(x)-abs(min(x))/10,max(y)+abs(max(y))/10,format_rho(rho));  xlim([min(x)-abs(min(x))/10,max(x)+abs(max(x))/10]); ylim([min(y)-abs(min(y))/10,max(y)+abs(max(y))/4]); 

subplot(4,6,16); 
plot(squeeze(mean(mean(mdb(:,5,alphabeta,other_orientations),1),3)),squeeze(mean(boldlow(other_orientations,:),2)),'d','LineWidth',1,'Color',[0.3,0.3,0.3]);  hold on ; 
plot(squeeze(mean(mean(mdb(:,5,alphabeta,oblique),1),3)),squeeze(mean(boldlow(oblique,:),2)),'d','LineWidth',2,'Color',oblique_color);
plot(squeeze(mean(mean(mdb(:,5,alphabeta,vertical),1),3)),squeeze(mean(boldlow(vertical,:),2)),'d','LineWidth',2,'Color',vertical_color); 
plot(squeeze(mean(mean(mdb(:,5,alphabeta,horizontal),1),3)),squeeze(mean(boldlow(horizontal,:),2)),'d','LineWidth',2,'Color',horizontal_color); 
x = double(squeeze(mean(mean(mdb(:,5,alphabeta,:),1),3))); y = double(squeeze(mean(boldlow,2))); 
B = [ones(size(x(:)))  x(:)]\y(:);
y_fit = [ones(size(x(:)))  x(:)]*B;
plot(x,y_fit,'k','LineWidth',3); xlabel('alpha/beta 5% contrast (dB)'); ylabel('BOLD 5% contrast (AU)'); title('5%contrast (alpha/beta vs BOLD)'); 
[rho,p] = corr(x,y); text(min(x)-abs(min(x))/10,min(y),format_rho(rho));  xlim([min(x)-abs(min(x))/10,max(x)+abs(max(x))/10]); ylim([min(y)-abs(min(y))/10,max(y)+abs(max(y))/4]); 
set(gcf,'Position',[100 100 1700 900])

subplot(4,6,17) ;
plot(1,'Color',oblique_color,'LineWidth',5); hold on; 
plot(1,'Color',vertical_color,'LineWidth',5);
plot(1,'Color',horizontal_color,'LineWidth',5);
plot(1,'Color',[0.3,0.3,0.3],'LineWidth',2);
legend({'oblique','vertical','horizontal','other orientations'});

figure,
subplot(1,2,1); 
plot(zscore(squeeze(mean(mean(allhigh(2:end,:,:),1),3)))); hold on ; 
plot(zscore(squeeze(mean(mean(mdb(:,1,gamma,:),1),3)))); 
subplot(1,2,2); 
plot(zscore(squeeze(mean(mean(alllow(2:end,:,:),1),3)))); hold on ; 
plot(zscore(squeeze(mean(mean(mdb(:,5,alphabeta,:),1),3)))); 


%}