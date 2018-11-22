clear all ; close all ; 
subs = {'lowcont_carmen','lowcont_lyndis','lowcont_mael','lowcont_lorne','lowcont_lyes','lowcont_greg','lowcont_angelina','lowcont_cesar','lowcont_charles','lowcont_alex','lowcont_janie','lowcont_russell','lowcont_anne'} ; 
comps = {[4,12,3],[6,40],[18,17],[8,6],[9,13],[5,6],[11,32],[38,4],[9,7],[22,16],[14,11],[27,17],[2,8]} ; 
stims = {'S  2','S  4','S 12','S 14','S 22','S 24','S 32','S 34','S 42','S 44'} ;

% methods figures: 
% use more frequency bands in the noise function 

cd e:/saved ; nrf = load('subersp') ; nrf = nrf.subersp; 
nrf_times = load('nrf_times'); nrf_times = nrf_times.nrf_times; 
nrf_freqs = load('nrf_freqs'); nrf_freqs = nrf_freqs.nrf_freqs; 

nrf_slowgamma = squeeze(mean(mean(nrf(:,nrf_freqs>25 & nrf_freqs<35,:),1),2)); 
nrf_fastgamma = squeeze(mean(mean(nrf(:,nrf_freqs>40 & nrf_freqs<60,:),1),2));
nrf_alphabeta = squeeze(mean(mean(nrf(:,nrf_freqs>10 & nrf_freqs<20,:),1),2));
nrf_low = squeeze(mean(mean(nrf(:,nrf_freqs>1 & nrf_freqs<6,:),1),2)); 

for sub=1:length(subs)
    
    cd(['E:/jly_orientation/',subs{sub}]) ;  

    merged = pop_loadset('merged.set') ;

    fullweights = load('fullweights'); fullweights = fullweights.fullweights; 
    eegweights = load('eegweights'); eegweights = eegweights.eegweights; 
      
    epmerged = merged; epmerged.data = eegweights{1}*eegweights{2}*merged.data; 
    
    
    
    
    %{
    clear mstersp; 
    for i=1:length(stims); disp(i); 
        allep = pop_epoch(epmerged,{stims{i}},[-2,20]); 
        for j=1:64
            for k=1:size(allep.data,3)
                [mstersp(i,j,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(j,:,k)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                        'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',60,'baseline',NaN,'verbose','off','timesout',200) ; 
            end
        end
    end
    bersp = mstersp - repmat(mean(mstersp(:,:,:,:,times<0),5),[1,1,1,1,200]);     
    resbesrp = zeros(size(bersp,2),size(bersp,1)*size(bersp,2),size(bersp,3),size(bersp,4)); 
    for i=1:10
       resbersp(:,((i-1)*size(bersp,3)+1):(i*size(bersp,3)),:,:) = squeeze(bersp(i,:,:,:,:));  
    end    
    [sv,si] = sort(squeeze(zscore(mean(mean(mean(abs(resbersp(:,:,35:end,:)),3),4),1))),'descend'); 
    goodtrials = si(10:end); 
    
    mbersp_slowgamma = squeeze(mean(bersp(:,:,:,nrf_freqs>25 & nrf_freqs<35,:),4)); 
    mbersp_fastgamma = squeeze(mean(bersp(:,:,:,nrf_freqs>40 & nrf_freqs<60,:),4)); 
    mbersp_alphabeta = squeeze(mean(bersp(:,:,:,nrf_freqs>10 & nrf_freqs<20,:),4)); 
    mbersp_low = squeeze(mean(bersp(:,:,:,nrf_freqs>1 & nrf_freqs<6,:),4)); 
    
    clear corrs_fastgamma corrs_slowgamma corrs_alphabeta corrs_low 
    
    for i=1:size(mbersp_fastgamma,1)
        for j=1:size(mbersp_fastgamma,2)
            for k=1:size(mbersp_fastgamma,3)
                corrs_fastgamma(i,j,k) = corr2(nrf_fastgamma,squeeze(mbersp_fastgamma(i,j,k,:))); 
                corrs_slowgamma(i,j,k) = corr2(nrf_slowgamma,squeeze(mbersp_slowgamma(i,j,k,:))); 
                corrs_alphabeta(i,j,k) = corr2(nrf_alphabeta,squeeze(mbersp_alphabeta(i,j,k,:))); 
                corrs_low(i,j,k) = corr2(nrf_low,squeeze(mbersp_low(i,j,k,:))); 

            end
        end
    end
    
    mcorrs = corrs_fastgamma + corrs_slowgamma + corrs_alphabeta/2 + corrs_low; 
    mcorrs(isnan(mcorrs)) = 0; 
    [sv,si] = sort(squeeze(mean(mean(mcorrs(1:8,:,:),1),3)),'descend'); 
    epoch_ica_nrfcorrs = si; save('epoch_ica_nrfcorrs','epoch_ica_nrfcorrs'); 
    figure,for i=1:64 ; subplot(5,13,i); imagesc(squeeze(mean(resbersp(si(i),:,:,:),2)),[-3,3]); title(i); axis xy; colormap jet; end; suptitle(subs{sub});
    
    subersp(sub,:,:) = squeeze(mean(mean(resbersp(comps{sub},goodtrials,:,:),1),2)); 
    %}
    
    
    epoch_ica_nrfcorrs = load('epoch_ica_nrfcorrs'); compcorrs = epoch_ica_nrfcorrs.epoch_ica_nrfcorrs; 
    clear mstersp; 
    for i=1:length(stims); disp(i); 
        allep = pop_epoch(epmerged,{stims{i}},[-2,20]); 
        for j=1:5
            for k=1:size(allep.data,3)
                [mstersp(i,j,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(compcorrs(j),:,k)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                        'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',80,'baseline',NaN,'verbose','off','timesout',1000) ; 
            end
        end
    end
    bersp = mstersp - repmat(mean(mstersp(:,:,:,:,times<0),5),[1,1,1,1,1000]);     
    mbersp = squeeze(mean(mean(mean(abs(bersp),2),4),5)); 
    [sv,si] = sort(mbersp,2,'descend'); 
    
    %figure,
    %for i=1:10 ; subplot(2,5,i) ; imagesc(squeeze(mean(mean(bersp(i,:,si(i,3:end),:,:),2),3)),[-2,2]) ; end
    
    for i=1:10; allbersp(sub,i,:,:) = squeeze(mean(mean(bersp(i,:,si(i,3:end),:,:),2),3)); end
    
    
        
end

time_inds = find(times > 0 & times < 18);    
gtime_inds = find(times>0.5 & times<18);    

oblique = [40:50,130:140];
vertical = 85:95; 
horizontal = [1:5,175:180];
cardinal = [85:95,1:5,175:180]; 

s2_angles = [90:-1:0,359:-1:45];
s4_angles = [45:1:359,0:1:90];

tbersp = allbersp(:,:,:,time_inds); 
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

tbersp = allbersp(:,:,:,gtime_inds); 
res_angles = res_angles(:,(length(time_inds)-length(gtime_inds)+1):end); 

tbersp_amp = zeros(size(tbersp,1),size(tbersp,2),size(tbersp,3),360);    
for i=0:359
    for j=1:10
        bersp_ind = find(i==res_angles(j,:)); 
        tbersp_amp(:,j,:,i+1) = squeeze(mean(tbersp(:,j,:,bersp_ind),4)); 
    end
end



%{
cd e:/saved; save('subersp','subersp');  
nrf_times = times; save('nrf_times','nrf_times'); 
nrf_freqs = freqs; save('nrf_freqs','nrf_freqs'); 
%}

clear sub_db
sub_db(:,1,:,:) = squeeze(mean(tbersp_amp(:,1:2,:,:),2)); 
sub_db(:,2,:,:) = squeeze(mean(tbersp_amp(:,3:4,:,:),2)); 
sub_db(:,3,:,:) = squeeze(mean(tbersp_amp(:,5:6,:,:),2)); 
sub_db(:,4,:,:) = squeeze(mean(tbersp_amp(:,7:8,:,:),2)); 
sub_db(:,5,:,:) = squeeze(mean(tbersp_amp(:,9:10,:,:),2)); 
mdb = (sub_db(:,:,:,1:180) + sub_db(:,:,:,181:end)) / 2;

% supplementary_01
subplot(3,5,1); imagesc(nrf_times,nrf_freqs,squeeze(mean(nrf)),[-1,1]) ; axis xy; title('neuronal response function (NRF)'); colormap jet; xlabel('time(s)'); ylabel('frequency(hz)'); 
subplot(3,5,2); 
shadedErrorBar(nrf_times,squeeze(mean(mean(nrf(:,nrf_freqs>10 & nrf_freqs<20,:),1),2)),squeeze(std(mean(nrf(:,nrf_freqs>10 & nrf_freqs<20,:),2),0,1))/sqrt(13),{'Color','b'}); hold on ; 
shadedErrorBar(nrf_times,squeeze(mean(mean(nrf(:,nrf_freqs>25 & nrf_freqs<35,:),1),2)),squeeze(std(mean(nrf(:,nrf_freqs>25 & nrf_freqs<35,:),2),0,1))/sqrt(13),{'Color','r'});
shadedErrorBar(nrf_times,squeeze(mean(mean(nrf(:,nrf_freqs>40 & nrf_freqs<60,:),1),2)),squeeze(std(mean(nrf(:,nrf_freqs>40 & nrf_freqs<60,:),2),0,1))/sqrt(13),{'Color','m'}); 
shadedErrorBar(nrf_times,squeeze(mean(mean(nrf(:,nrf_freqs>1 & nrf_freqs<6,:),1),2)),squeeze(std(mean(nrf(:,nrf_freqs>1 & nrf_freqs<6,:),2),0,1))/sqrt(13),{'Color','c'}); xlim([-2,20]); xlabel('time(s)'); ylabel('dB'); 
vline([0,18]); hline(0,'k'); xlim([-2,20]); ylabel('dB'); xlabel('time(s)'); title('frequency specific NRF'); 
subplot(3,5,3); imagesc([-1,1]) ; h = colorbar ; title(h,'dB'); 
subplot(3,5,4); plot(1,'c'); hold on ; plot(1,'b'); plot(1,'r'); plot(1,'m'); legend({'1-6hz','10-20hz','25-35hz','40-60hz'});

subplot(3,5,6) ; imagesc(nrf_times,nrf_freqs,squeeze(mean(allbersp(:,1,:,:),1)),[-3,3]); axis xy; colormap jet; xlabel('time(s)'); ylabel('frequency(hz)'); title('clockwise');
subplot(3,5,7) ; imagesc(nrf_times,nrf_freqs,squeeze(mean(allbersp(:,2,:,:),1)),[-3,3]); axis xy; colormap jet; xlabel('time(s)'); ylabel('frequency(hz)'); title('counter-clockwise'); 
subplot(3,5,8) ; imagesc(1:360,freqs,squeeze(mean(sub_db(:,1,:,:),1)),[-3,3]);  axis xy; xlabel('orientation (deg\circ)'); title('1-360\circ'); set(gca,'XTick',[0,90,180,270,360]); 
subplot(3,5,9); imagesc(1:180,freqs,squeeze(mean(mdb(:,1,:,:),1)),[-3,3]); axis xy; xlabel('orientation (deg\circ)'); title('1-180\circ'); set(gca,'XTick',[45,90,135,180]); 
subplot(3,5,10) ; imagesc([-3,3]); h = colorbar ; title(h,'dB'); 
set(gcf,'Position',[100 100 1700 900])


[xg,yg] = meshgrid([-100:.25:100,100:.25:100]); 
circmask = sqrt(xg.^2 + yg.^2) <100; 
redxhair = sqrt(xg.^2 + yg.^2) <4; 
stim_contrasts = [1,.5,.25,.15,.05];
for i=1:length(stim_contrasts)
    contrast_i = sin(xg); 
    contrast_i = contrast_i*122.5*stim_contrasts(i) + 122.5;
    x = uint8(contrast_i.*circmask); x(circmask==0) = 122.5; 
    subplot(2,3,i); imshow(uint8(x)); 
    title([num2str(stim_contrasts(i)*100),'%']);
end

rotangles = [-90,-45,0,45,90,135,180,225,270];

[xg,yg] = meshgrid([-100:.25:100,100:.25:100]); 
circmask = sqrt(xg.^2 + yg.^2) < 50; 
contrast_i = sin(xg); 
contrast_i = contrast_i*122.5*stim_contrasts(1) + 122.5;
for i=1:9 
    subplot(3,3,i) ; 
    x = uint8(imrotate(contrast_i.*circmask,rotangles(i),'crop')); x(circmask==0) = 255; 
    imshow(x(end/4:end-end/4,end/4:end-end/4)); title([num2str(rotangles(i)+90),'\circ']); 
end


cd e:/saved; 
allhigh = load('allhigh'); allhigh = allhigh.allhigh; 
alllow = load('alllow'); alllow = alllow.alllow; 

contrasts = {'100%','30%','15%','5%','1%'};
for i=1:5 ; subplot(3,7,i) ; 
    imagesc(1:360,freqs,squeeze(mean(mdb(:,i,:,:),1))) ; axis xy ; colormap jet; if i==1; xlabel('bar orientation (deg\circ)'); ylabel('frequency(hz)');  end; title(contrasts{i});
    h = colorbar ; title(h,'dB'); 
end 

sfreqs = {20:35,14:16,7:12,4:6};
freqtitles = {'mid gamma (40-70Hz)','low gamma (28-32Hz)','beta (15-25Hz)','alpha (8-12Hz)'};
freq_colors = {[0,0,0.5],[0,0.5,0],[0.5,0,0],[0.25,0.15,0.15]};
contrast_labs = {'100','30','15','5','1'};

% row 2 - show all amplitude cardinal vs oblique
for i=1:length(freqtitles)
    subplot(3,7,i+7); 
    for j=1:5
    card_amp = squeeze(mean(mean(mdb(:,j,sfreqs{i},cardinal),3),4)); 
    obl_amp = squeeze(mean(mean(mdb(:,j,sfreqs{i},oblique),3),4)); 
    
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
   [h,p,ci,stats] = ttest(squeeze(mean(mean(mdb(:,i,j,cardinal),4),3)),squeeze(mean(mean(mdb(:,i,j,oblique),4),3))); 
   allts(i,j) = stats.tstat; allps(i,j) = p; 
    end
end

% PEAK FREQUENCY ESTIMATION:
% zoom in on gamma, plot all subjects with peak frequency (low and mid)
mspec = squeeze(mean(mean(mean(mdb,1),2),4)); 
contrast_colors = {[0.1,0.1,0],[0.2,0.2,0],[0.3,0.3,0],[0.4,0.4,0],[0.5,0.5,0]}; 
oblique_color = 'r'; vertical_color = 'c'; horizontal_color = 'g'; cardinal_color = 'm'; 

ideal_mid = fspecial('gaussian',[1,15],3); 
ideal_low = fspecial('gaussian',[1,10],2); 

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
   shadedErrorBar(freqs(12:35),squeeze(mean(mean(mdb(:,i,12:35,:),1),4)),squeeze(std(mean(mdb(:,i,12:35,:),4),0,1))/sqrt(13),{'Color',contrast_colors{i}}); hold on ;
   xlabel('frequency(hz)'); ylabel('dB'); xlim([23,70]);
    
end
title('mean gamma (all orientations)'); 
subplot(4,3,3) ; for i=1:5 ; plot(1,'Color',contrast_colors{i},'LineWidth',5) ; hold on ; end ; legend({'100% contrast','50% contrast','25% contrast','15% contrast','5% contrast'});

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

subplot(4,6,16) ; plot(1,oblique_color,'LineWidth',5) ; hold on ; plot(1,cardinal_color,'LineWidth',5); legend({'oblique','cardinal'});
subplot(4,6,17) ; plot(1,oblique_color,'LineWidth',5) ; hold on ; plot(1,vertical_color,'LineWidth',5); legend({'oblique','vertical'});
subplot(4,6,18) ; plot(1,vertical_color,'LineWidth',5) ; hold on ; plot(1,horizontal_color,'LineWidth',5); legend({'vertical','horizontal'});


set(gcf,'Position',[100 100 1700 900])



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

%%% figure 2: eeg orientation tuning
mdb = (sub_db(:,:,:,1:180) + sub_db(:,:,:,181:end)) / 2;


% sort the orientations
for i=1:13
    for j=1:5
        [~,sort_gammadb(i,j,:)] = sort(squeeze(mean(mdb(i,j,14:35,:),3)),'descend'); 
        [~,sort_alphadb(i,j,:)] = sort(squeeze(mean(mdb(i,j,6:8,:),3)),'descend'); 
    end
end

orientation_values = {oblique,cardinal,vertical,horizontal};

for i=1:13
    for j=1:5
        gamma_ij = squeeze(sort_gammadb(i,j,:)); 
        alpha_ij = squeeze(sort_alphadb(i,j,:)); 
        top_gamma_oblique(i,j) = length(intersect(oblique,gamma_ij(1:40)))/40; 
        top_gamma_cardinal(i,j) = length(intersect(cardinal,gamma_ij(1:40)))/40; 
        top_alpha_vertical(i,j) = length(intersect(vertical,alpha_ij(1:40)))/40; 
        top_alpha_horizontal(i,j) = length(intersect(horizontal,alpha_ij(1:40)))/40; 
    end
end



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
%{
for i=1:13  
   zgl(i,:) = (zscore(squeeze(mean(mdb(i,5,gamma,:),3))));   
   zgh(i,:) = (zscore(squeeze(mean(mdb(i,1,gamma,:),3))));   
   zal(i,:) = (zscore(squeeze(mean(mdb(i,5,alphabeta,:),3))));   
   zah(i,:) = (zscore(squeeze(mean(mdb(i,1,alphabeta,:),3))));   
end
subplot(4,6,8);
shadedErrorBar([],squeeze(mean(zgl,1)),squeeze(std(zgl,0,1))/sqrt(13),{'Color',contrast_colors{5}}); hold on; 
shadedErrorBar([],squeeze(mean(zgh,1)),squeeze(std(zgh,0,1))/sqrt(13),{'Color',contrast_colors{1}}); ylabel('z-score'); xlabel('orientation(deg\circ)'); 
subplot(4,6,7); 
shadedErrorBar([],squeeze(mean(zal,1)),squeeze(std(zal,0,1))/sqrt(13),{'Color',contrast_colors{5}}); hold on; 
shadedErrorBar([],squeeze(mean(zah,1)),squeeze(std(zah,0,1))/sqrt(13),{'Color',contrast_colors{1}}); ylabel('z-score'); xlabel('orientation(deg\circ)'); 
%}


subplot(4,6,7); 
shadedErrorBar(1:180,squeeze(mean(mean(mdb(:,1,alphabeta,:),1),3)),squeeze(std(mean(mdb(:,1,alphabeta,:),3),0,1))/sqrt(13),{'Color',contrast_colors{1}}); hold on ; 
shadedErrorBar(1:180,squeeze(mean(mean(mdb(:,5,alphabeta,:),1),3)),squeeze(std(mean(mdb(:,5,alphabeta,:),3),0,1))/sqrt(13),{'Color',contrast_colors{5}}); hold on ; 
set(gca,'XTick',[45,90,135]) ; xlabel('orientation (deg\circ)'); ylabel('dB'); grid on; title('alpha/beta orientation tuning, *p<0.05'); ylim([-3.5,-1.4]);
for i=1:180 ; [~,p,~,~] = ttest(squeeze(mean(mdb(:,1,alphabeta,i),3)),squeeze(mean(mdb(:,5,alphabeta,i),3))); if p<0.05 ; text(i,-1.5,'*') ; end; end ; xlim([0,180]);

subplot(4,6,8); 
shadedErrorBar(1:180,squeeze(mean(mean(mdb(:,1,gamma,:),1),3)),squeeze(std(mean(mdb(:,1,gamma,:),3),0,1))/sqrt(13),{'Color',contrast_colors{1}}); hold on ; 
shadedErrorBar(1:180,squeeze(mean(mean(mdb(:,5,gamma,:),1),3)),squeeze(std(mean(mdb(:,5,gamma,:),3),0,1))/sqrt(13),{'Color',contrast_colors{5}}); hold on ; 
set(gca,'XTick',[45,90,135]) ; xlabel('orientation (deg\circ)'); ylabel('dB'); grid on; title('gamma orientation tuning, *p<0.05'); ylim([-.5,1.9]);
for i=1:180 ; [~,p,~,~] = ttest(squeeze(mean(mdb(:,1,gamma,i),3)),squeeze(mean(mdb(:,5,gamma,i),3))); if p<0.05 ; text(i,1.8,'*') ; end; end;xlim([0,180]);


subplot(4,12,17) ; plot(1,'Color',contrast_colors{1},'LineWidth',5) ; hold on ; plot(1,'Color',contrast_colors{5},'LineWidth',5) ; legend({'100% contrast','5% contrast'});



subplot(4,6,10); 
for i=1:5
   bar(i*3-1,squeeze(mean(top_gamma_oblique(:,i))),'r'); hold on;    
   bar(i*3,squeeze(mean(top_gamma_cardinal(:,i))),'m');  
   errorbar((1:3:15)+1,squeeze(mean(top_gamma_oblique,1)),squeeze(std(top_gamma_oblique,0,1))/sqrt(13),'.k'); 
   errorbar((2:3:15)+1,squeeze(mean(top_gamma_cardinal,1)),squeeze(std(top_gamma_cardinal,0,1))/sqrt(13),'.k'); 
   maxvals(i) = max(squeeze(mean(top_gamma_oblique(:,i),1)) + squeeze(std(top_gamma_oblique(:,i),0,1))/sqrt(13)); 
end
ylabel('proportion of top 20%'); xlabel('contrast(%)');
for i=1:5
    [h,ps,ci,stats] = ttest(top_gamma_oblique(:,i),top_gamma_cardinal(:,i));  
    if ps<0.05 ;text(i*3-.5,maxvals(i),'*'); end    
end
title('*p<0.05'); set(gca,'XTick',2.5:3:15,'XTickLabel',[100,50,25,15,5]); 

subplot(4,6,11); 
for i=1:5
   bar(i*3-1,squeeze(mean(top_alpha_vertical(:,i))),'c'); hold on;    
   bar(i*3,squeeze(mean(top_alpha_horizontal(:,i))),'g');   
   errorbar((1:3:15)+1,squeeze(mean(top_alpha_vertical,1)),squeeze(std(top_alpha_vertical,0,1))/sqrt(13),'.k'); 
   errorbar((2:3:15)+1,squeeze(mean(top_alpha_horizontal,1)),squeeze(std(top_alpha_horizontal,0,1))/sqrt(13),'.k');    
   maxvals(i) = max(squeeze(mean(top_alpha_horizontal(:,i),1)) + squeeze(std(top_alpha_horizontal(:,i),0,1))/sqrt(13)); 
end
ylabel('proportion of top 20%'); xlabel('contrast(%)');
for i=1:5
    [h,ps,ci,stats] = ttest(top_alpha_vertical(:,i),top_alpha_horizontal(:,i));  
    if ps<0.05 ;text(i*3-.5,maxvals(i),'*'); end    
end
title('*p<0.05'); set(gca,'XTick',2.5:3:15,'XTickLabel',[100,50,25,15,5]); 


%{
clear mdb_corrs; 
for i=1:5
   mdb_corrs(i,:,:) = corr(squeeze(mean(mdb(:,i,:,:),1))',squeeze(mean(mdb(:,1,:,:),1))');  
   for j=1:size(mdb_corrs,2)
       for k=1:size(mdb_corrs,3)
         if j<k ; mdb_corrs(i,j,k) = 0; end
       end
   end
end
titles = {'100%','50%','25%','15%','5%'};
for i=1:5 ;subplot(4,12,18+i) ; imagesc(freqs,freqs,squeeze(mdb_corrs(i,:,:)),[-.8,.8]); colormap jet; title(['100% vs ',titles{i}]); if i==1; xlabel('frequency(hz)'); ylabel(['frequency(hz)']); end; end
subplot(4,12,24) ; imagesc([-.8,.8]) ;h=colorbar; title(h,'rho'); set(gca,'XTick',[],'YTick',[]); 
%}

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

subplot(4,6,16) ; plot(1,oblique_color,'LineWidth',5); hold on ; plot(1,cardinal_color,'LineWidth',5); legend({'oblique','cardinal'});
subplot(4,6,17) ; plot(1,oblique_color,'LineWidth',5); hold on ; plot(1,vertical_color,'LineWidth',5); legend({'oblique','vertical'});
subplot(4,6,18) ; plot(1,vertical_color,'LineWidth',5); hold on ; plot(1,horizontal_color,'LineWidth',5); legend({'vertical','horizontal'});
set(gcf,'Position',[100 100 1700 900])



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

%{
subplot(4,3,1) ; 
imagesc(freqs,1:5,highcorrs,[-1,1]) ; colormap jet; xlabel('frequency(hz)') ; ylabel('EEG contrast(%)') ;set(gca,'YTick',1:5,'YTickLabel',[100,50,25,15,5]); hline([0.5:5.5],'k'); 
title('EEG (all contrast) vs BOLD (100% contrast)'); 
subplot(4,3,2); 
imagesc(freqs,1:5,lowcorrs,[-1,1]) ; colormap jet; xlabel('frequency(hz)') ; ylabel('EEG contrast(%)') ;set(gca,'YTick',1:5,'YTickLabel',[100,50,25,15,5]); hline([0.5:5.5],'k'); 
title('EEG (all contrast) vs BOLD (5% contrast)'); 
subplot(4,3,3) ; imagesc([-1,1]); h = colorbar; title(h,'rho'); 
%}

subplot(4,3,1); 
plot(zscore(squeeze(mean(mean(allhigh(2:end,:,:),1),3)))); hold on ; 
plot(zscore(squeeze(mean(mean(mdb(:,1,gamma,:),1),3)))); 
legend({'BOLD','EEG gamma'});xlabel('orientation(deg\circ)'); ylabel('tuning z-score'); xlim([0,180]) ; set(gca,'XTick',[45,90,135]); 
title('100% contrast'); 
subplot(4,3,2); 
plot(zscore(squeeze(mean(mean(alllow(2:end,:,:),1),3)))); hold on ; 
plot(zscore(squeeze(mean(mean(mdb(:,5,alphabeta,:),1),3)))); 
legend({'BOLD','EEG alpha/beta'});xlabel('orientation(deg\circ)'); ylabel('tuning z-score'); xlim([0,180]) ; set(gca,'XTick',[45,90,135]); 
title('5% contrast'); 

orientations = zeros(1,180) ; orientations([oblique,vertical,horizontal]) = 1; other_orientations = find(orientations==0); 

subplot(4,6,7); 
plot(squeeze(mean(mean(mdb(:,1,gamma,other_orientations),1),3)),squeeze(mean(boldhigh(other_orientations,:),2)),'d','LineWidth',1,'Color',[0.3,0.3,0.3]);  hold on ; 
plot(squeeze(mean(mean(mdb(:,1,gamma,oblique),1),3)),squeeze(mean(boldhigh(oblique,:),2)),'d','LineWidth',2,'Color',oblique_color);
plot(squeeze(mean(mean(mdb(:,1,gamma,vertical),1),3)),squeeze(mean(boldhigh(vertical,:),2)),'d','LineWidth',2,'Color',vertical_color); 
plot(squeeze(mean(mean(mdb(:,1,gamma,horizontal),1),3)),squeeze(mean(boldhigh(horizontal,:),2)),'d','LineWidth',2,'Color',horizontal_color); 
x = double(squeeze(mean(mean(mdb(:,1,gamma,:),1),3))); y = double(squeeze(mean(boldhigh,2))); 
B = [ones(size(x(:)))  x(:)]\y(:);
y_fit = [ones(size(x(:)))  x(:)]*B;
plot(x,y_fit,'k','LineWidth',3); xlabel('gamma 100% contrast (dB)'); ylabel('BOLD 100% contrast (AU)'); title('100%contrast (gamma vs BOLD)'); 
[rho,p] = corr(x,y); text(min(x)-abs(min(x))/10,max(y)+abs(max(y))/10,format_rho(rho)); xlim([min(x)-abs(min(x))/10,max(x)+abs(max(x))/10]); ylim([min(y)-abs(min(y))/10,max(y)+abs(max(y))/4]); 

subplot(4,6,8); 
plot(squeeze(mean(mean(mdb(:,1,alphabeta,other_orientations),1),3)),squeeze(mean(boldhigh(other_orientations,:),2)),'d','LineWidth',1,'Color',[0.3,0.3,0.3]);  hold on ; 
plot(squeeze(mean(mean(mdb(:,1,alphabeta,oblique),1),3)),squeeze(mean(boldhigh(oblique,:),2)),'d','LineWidth',2,'Color',oblique_color);
plot(squeeze(mean(mean(mdb(:,1,alphabeta,vertical),1),3)),squeeze(mean(boldhigh(vertical,:),2)),'d','LineWidth',2,'Color',vertical_color); 
plot(squeeze(mean(mean(mdb(:,1,alphabeta,horizontal),1),3)),squeeze(mean(boldhigh(horizontal,:),2)),'d','LineWidth',2,'Color',horizontal_color); 
x = double(squeeze(mean(mean(mdb(:,1,alphabeta,:),1),3))); y = double(squeeze(mean(boldhigh,2))); 
B = [ones(size(x(:)))  x(:)]\y(:);
y_fit = [ones(size(x(:)))  x(:)]*B;
plot(x,y_fit,'k','LineWidth',3); xlabel('alpha/beta 100% contrast (dB)'); ylabel('BOLD 100% contrast (AU)'); title('100%contrast (alpha/beta vs BOLD)'); 
[rho,p] = corr(x,y); text(min(x)-abs(min(x))/10,max(y)+abs(max(y))/10,format_rho(rho));   xlim([min(x)-abs(min(x))/10,max(x)+abs(max(x))/10]); ylim([min(y)-abs(min(y))/10,max(y)+abs(max(y))/4]); 

subplot(4,6,9); 
plot(squeeze(mean(mean(mdb(:,5,gamma,other_orientations),1),3)),squeeze(mean(boldlow(other_orientations,:),2)),'d','LineWidth',1,'Color',[0.3,0.3,0.3]);  hold on ; 
plot(squeeze(mean(mean(mdb(:,5,gamma,oblique),1),3)),squeeze(mean(boldlow(oblique,:),2)),'d','LineWidth',2,'Color',oblique_color);
plot(squeeze(mean(mean(mdb(:,5,gamma,vertical),1),3)),squeeze(mean(boldlow(vertical,:),2)),'d','LineWidth',2,'Color',vertical_color); 
plot(squeeze(mean(mean(mdb(:,5,gamma,horizontal),1),3)),squeeze(mean(boldlow(horizontal,:),2)),'d','LineWidth',2,'Color',horizontal_color); 
x = double(squeeze(mean(mean(mdb(:,5,gamma,:),1),3))); y = double(squeeze(mean(boldlow,2))); 
B = [ones(size(x(:)))  x(:)]\y(:);
y_fit = [ones(size(x(:)))  x(:)]*B;
plot(x,y_fit,'k','LineWidth',3); xlabel('gamma 5% contrast (dB)'); ylabel('BOLD 5% contrast (AU)'); title('5%contrast (gamma vs BOLD)'); 
[rho,p] = corr(x,y); text(min(x)-abs(min(x))/10,max(y)+abs(max(y))/10,format_rho(rho));  xlim([min(x)-abs(min(x))/10,max(x)+abs(max(x))/10]); ylim([min(y)-abs(min(y))/10,max(y)+abs(max(y))/4]); 

subplot(4,6,10); 
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

subplot(4,6,12) ;
plot(1,'Color',oblique_color,'LineWidth',5); hold on; 
plot(1,'Color',vertical_color,'LineWidth',5);
plot(1,'Color',horizontal_color,'LineWidth',5);
plot(1,'Color',[0.3,0.3,0.3],'LineWidth',2);
legend({'oblique','vertical','horizontal','other orientations'});

subplot(4,3,7) ; 
plot(freqs,squeeze(highcorrs(1,:)),'Color',contrast_colors{1},'LineWidth',2); hold on ; hline(0,'k'); 
plot(freqs,squeeze(lowcorrs(5,:)),'Color',contrast_colors{5},'LineWidth',2); hold on ; hline(0,'k'); 
xlim([0,120]); ylabel('rho'); xlabel('frequency(hz)'); title('EEG-BOLD tuning correlation'); 

subplot(4,3,8) ; 
plot(freqs,squeeze(highcorrs(1,:)).^2,'Color',contrast_colors{1},'LineWidth',2); hold on ; hline(0,'k'); 
plot(freqs,squeeze(lowcorrs(5,:)).^2,'Color',contrast_colors{5},'LineWidth',2); hold on ; hline(0,'k'); 
xlim([0,120]); ylabel('rho^{2}'); xlabel('frequency(hz)'); title('BOLD tuning variance explained'); 

subplot(4,3,9); plot(1,'Color',contrast_colors{1},'LineWidth',5); hold on; plot(1,'Color',contrast_colors{5},'LineWidth',5); legend({'100% contrast','5% contrast'});

msubalpha_high = squeeze(mean(mean(mdb(:,1,freqs>10 & freqs<18,:),1),3)); 
msubgamma_high = squeeze(mean(mean(mdb(:,1,freqs>28 & freqs<60,:),1),3)); 
msubalpha_low = squeeze(mean(mean(mdb(:,5,freqs>10 & freqs<18,:),1),3)); 
msubgamma_low = squeeze(mean(mean(mdb(:,5,freqs>28 & freqs<60,:),1),3)); 
mboldhigh = squeeze(mean(mean(allhigh,1),3)); 
mboldlow = squeeze(mean(mean(alllow,1),3)); 

bothgamma = [msubgamma_low;msubgamma_high];
bothalpha = [msubalpha_low;msubalpha_high];
bothbold = [mboldlow,mboldhigh];
botheeg = [bothgamma,bothalpha];

[beta,sigma,resid] = mvregress(botheeg,bothbold'); 








