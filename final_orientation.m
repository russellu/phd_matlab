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
    %fullweights = load('fullweights'); fullweights = fullweights.fullweights; 
    eegweights = load('eegweights'); eegweights = eegweights.eegweights; 

    
     %merged = pop_loadset('merged2.set') ;
     %eegweights = load('eegcomps2'); eegweights = eegweights.eegcomps2; 
    
    
    % get the scalp maps 
    
    clocs{sub} = merged.chanlocs; 
    winv = pinv(eegweights{1}*eegweights{2});  
    epoch_ica_nrfcorrs = load('epoch_ica_nrfcorrs'); compcorrs = epoch_ica_nrfcorrs.epoch_ica_nrfcorrs; 
    %compcorrs = load('newcomps_si') ; compcorrs = compcorrs.newcomps_si; 
    epmerged = merged; epmerged.data = eegweights{1}*eegweights{2}*merged.data; 
    epmerged.data(compcorrs(6:end),:) = 0; 
    epmerged.data = winv*epmerged.data; 
    clear ctrigs
    ctrigs{1} = {'S  2','S  4'}; ctrigs{2} = {'S 12','S 14'}; ctrigs{3} = {'S 22','S 24'}; ctrigs{4} = {'S 32','S 34'}; ctrigs{5} = {'S 42','S 44'}; 
    gamma_hz = 35:70; alpha_hz = 9:15; clear filtgamma filtalpha
    for i=1:length(ctrigs)
        epsi = pop_epoch(epmerged,ctrigs{i},[-2,20]); 
        for j=1:size(epsi.data,3)
            filtgamma(:,:,i,j) = (eegfiltfft(epsi.data(:,:,j),epsi.srate,40,65)); 
            filtalpha(:,:,i,j) = (eegfiltfft(epsi.data(:,:,j),epsi.srate,9,15)); 
        end
    end
    clear stds 
    for i=1:5
        for j=1:size(filtgamma,4)
            stds(i,j) = mean(std(filtgamma(:,:,i,j),0,2),1); 
        end
    end
    [sv,si] = sort(stds,2,'descend'); 
    
    filtgamma = abs(filtgamma);
    filtalpha = abs(filtalpha); 
    
    clear bfiltgamma bfiltalpha; 
    for i=1:5
    bfiltgamma(:,i,:) = squeeze(mean(filtgamma(:,epsi.times>0 & epsi.times<18000,i,si(i,5:end)),2) - mean(filtgamma(:,(epsi.times<0 & epsi.times>-1500),i,si(i,5:end)),2)); 
    bfiltalpha(:,i,:) = squeeze(mean(filtalpha(:,epsi.times>0 & epsi.times<18000,i,si(i,5:end)),2) - mean(filtalpha(:,(epsi.times<0 & epsi.times>-1500),i,si(i,5:end)),2)); 
    end
    
    gtopo(sub,:,:) = mean(bfiltgamma,3);
    atopo(sub,:,:) = mean(bfiltalpha,3);

    epmerged = merged; epmerged.data = eegweights{1}*eegweights{2}*merged.data; 
    epoch_ica_nrfcorrs = load('epoch_ica_nrfcorrs'); compcorrs = epoch_ica_nrfcorrs.epoch_ica_nrfcorrs; 
    %compcorrs = load('newcomps_si') ; compcorrs = compcorrs.newcomps_si; 
    
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
    for i=1:10; allbersp(sub,i,:,:) = squeeze(mean(mean(bersp(i,:,si(i,3:end),:,:),2),3)); end
    
         
end

clear gammagrid ; clear alphagrid; 
for i=1:13 
    for s=1:5
    subplot(5,13,i) ;     
    if i==7 || i==2
        [~,gammagrid(i,s,:,:),~] = topoplot(gtopo(i,:,s),clocs{1},'gridscale',200); 
        [~,alphagrid(i,s,:,:),~] = topoplot(atopo(i,:,s),clocs{1},'gridscale',200); 

    else
        [~,gammagrid(i,s,:,:),~] = topoplot(gtopo(i,:,s),clocs{i},'gridscale',200);        
        [~,alphagrid(i,s,:,:),~] = topoplot(atopo(i,:,s),clocs{i},'gridscale',200);        
    end   
    end
end


figure,
clabs = {'100%','50%','25%','15%','5%'};
cmap = jet(256); 
for i=1:5
nangamma_i = squeeze(mean(gammagrid(:,i,:,:),1)); 
gamma_i = uint8((nangamma_i)*1700); 
rgbgamma = zeros(size(gamma_i,1),size(gamma_i,2),3); 
nanlpha_i = squeeze(mean(alphagrid(:,i,:,:),1)); 
alpha_i = 255-uint8((-nanlpha_i)*127.5); 
rgbalpha = zeros(size(alpha_i,1),size(alpha_i,2),3); 
for x=1:size(nangamma_i,1)
    for y=1:size(nangamma_i,2)
        if ~ isnan(nangamma_i(x,y))
           rgbgamma(x,y,:) = cmap(gamma_i(x,y)+1,:); 
           rgbalpha(x,y,:) = cmap(alpha_i(x,y)+1,:); 
        else
           rgbgamma(x,y,:) = [255,255,255];  
           rgbalpha(x,y,:) = [255,255,255]; 
        end
        
    end
end
subplot(2,12,i);
imshow(rgbgamma) ; axis xy ; hold on; contour(gamma_i,6,'k'); title(clabs{i}); 
subplot(2,12,i+12);
imshow(rgbalpha) ; axis xy ; hold on; contour(alpha_i,6,'k'); 
end

subplot(5,5,5); imagesc([0,.2]) ; h = colorbar; title(h,'\muV (40-60HZ)'); colormap jet 
subplot(5,5,10); imagesc([-2,0]) ; h = colorbar; title(h,'\muV (10-16Hz)'); 
figure, topoplot([],merged.chanlocs,'style','blank','electrodes','labelpoint','chaninfo',merged.chaninfo);



time_inds = find(times > 0 & times < 18);    
gtime_inds = find(times>0.5 & times<18);    

oblique = [35:55,125:145];
vertical = 80:100; 
horizontal = [1:10,170:180];
cardinal = [80:100,1:10,170:180]; 

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

clear sub_db
sub_db(:,1,:,:) = squeeze(mean(tbersp_amp(:,1:2,:,:),2)); 
sub_db(:,2,:,:) = squeeze(mean(tbersp_amp(:,3:4,:,:),2)); 
sub_db(:,3,:,:) = squeeze(mean(tbersp_amp(:,5:6,:,:),2)); 
sub_db(:,4,:,:) = squeeze(mean(tbersp_amp(:,7:8,:,:),2)); 
sub_db(:,5,:,:) = squeeze(mean(tbersp_amp(:,9:10,:,:),2)); 
mdb = (sub_db(:,:,:,1:180) + sub_db(:,:,:,181:end)) / 2;

% PEAK FREQUENCY ESTIMATION:
% zoom in on gamma, plot all subjects with peak frequency (low and mid)
mspec = squeeze(mean(mean(mean(mdb,1),2),4)); 
CONTRAST_COLORS = {[194,82,60]/255,[230,142,28]/255,[123,237,0]/255,[30,144,148]/255,[11,44,122]/255}; 
CONTRAST_VALUES = [100,50,25,15,5]; 


ideal_mid = fspecial('gaussian',[1,15],3); 
ideal_low = fspecial('gaussian',[1,10],2); 
ideal_alpha = -ideal_low; 

orientation_inds = {1:180,oblique,cardinal,vertical,horizontal};

HZ_MIDGAMMA = find(freqs>40 & freqs<60); HZ_LOWGAMMA = find(freqs>25 & freqs<35) ; HZ_ALPHABETA = find(freqs>10 & freqs<16); 
HZ_BOTHGAMMA = find(freqs>25 & freqs<60); 
HZ_ARRAY = {HZ_MIDGAMMA,HZ_LOWGAMMA,HZ_ALPHABETA};


clear xcorrs lowxcorrs
for i=1:size(mdb,1)
    for j=1:size(mdb,2)
        for k=1:length(orientation_inds)
            xcorrs(i,j,k,:) = xcorr(squeeze(mean(mdb(i,j,15:40,orientation_inds{k}),4)),ideal_mid,25); 
            lowxcorrs(i,j,k,:) = xcorr(squeeze(mean(mdb(i,j,10:20,orientation_inds{k}),4)),ideal_low,12); 
            alphacorrs(i,j,k,:) = xcorr(squeeze(mean(mdb(i,j,2:12,orientation_inds{k}),4)),ideal_alpha,12); 
        end
    end
end

clear peak_midgamma peak_lowgamma; 
for i=1:size(mdb,1)
    for j=1:size(mdb,2)
        for k=1:length(orientation_inds)
            peak_midgamma(i,j,k) = find(max(squeeze(xcorrs(i,j,k,:)))==squeeze(xcorrs(i,j,k,:)),1); 
            peak_lowgamma(i,j,k) = find(max(squeeze(lowxcorrs(i,j,k,:)))==squeeze(lowxcorrs(i,j,k,:)),1); 
            peak_alpha(i,j,k) = find(max(squeeze(alphacorrs(i,j,k,:)))==squeeze(alphacorrs(i,j,k,:)),1);
        end
    end
end
peak_midgamma = (peak_midgamma - 27)*2 + 45;
peak_lowgamma = (peak_lowgamma - 11)*2 + 25; 
peak_alpha = (peak_alpha - 7)*2 + 1; 

figure,
for i=1:5 
   subplot(5,7,i) ; imagesc(1:180,freqs,squeeze(mean(mdb(:,i,:,:),1)),[-3.5,1.5]) ; axis xy ; colormap jet;  
   set(gca,'XTick',[45,90,135,180]); xlabel('orientation(deg\circ)'); ylabel('frequency(hz)'); title([num2str(CONTRAST_VALUES(i)), '% contrast']); 
end
subplot(5,7,7) ; imagesc([-3.5,1.5]) ; h = colorbar; title(h,'dB'); 

figure,
% mid-gamma amplitude
subplot(6,8,1);
for i=1:5
    h = bar(i*2,squeeze(mean(mean(mean(mdb(:,i,HZ_MIDGAMMA,:),1),3),4))); hold on ; 
    h.FaceColor = CONTRAST_COLORS{i}; 
    errorbar(i*2,squeeze(mean(mean(mean(mdb(:,i,HZ_MIDGAMMA,:),1),3),4)),squeeze(std(mean(mean(mdb(:,i,HZ_MIDGAMMA,:),4),3),0,1))/sqrt(13),'k.'); 
    [h,p,ci,stats] = ttest(squeeze(mean(mean(mdb(:,1,HZ_MIDGAMMA,:),3),4)),squeeze(mean(mean(mdb(:,i,HZ_MIDGAMMA,:),3),4))); 
    if p < 0.05 ; text(i*2,1,'*') ; end
end
set(gca,'XTick',2:2:10,'XTickLabel',CONTRAST_VALUES); xlabel('contrast(%)'); ylabel('dB'); title('mid gamma, *p<0.05'); xlim([0,13]); ylim([-.5,1.6]);
% low-gamma amplitude
subplot(6,8,2);
for i=1:5
    h = bar(i*2,squeeze(mean(mean(mean(mdb(:,i,HZ_LOWGAMMA,:),1),3),4))); hold on ; 
    h.FaceColor = CONTRAST_COLORS{i}; 
    errorbar(i*2,squeeze(mean(mean(mean(mdb(:,i,HZ_LOWGAMMA,:),1),3),4)),squeeze(std(mean(mean(mdb(:,i,HZ_LOWGAMMA,:),4),3),0,1))/sqrt(13),'k.'); 
    [h,p,ci,stats] = ttest(squeeze(mean(mean(mdb(:,1,HZ_LOWGAMMA,:),3),4)),squeeze(mean(mean(mdb(:,i,HZ_LOWGAMMA,:),3),4))); 
    if p < 0.05 ; text(i*2,1,'*') ; end
end
set(gca,'XTick',2:2:10,'XTickLabel',CONTRAST_VALUES); xlabel('contrast(%)'); ylabel('dB'); title('low gamma, *p<0.05'); xlim([0,13]);ylim([-.5,1.6]);
% alpha/beta amplitude
subplot(6,8,3);
for i=1:5
    h = bar(i*2,squeeze(mean(mean(mean(mdb(:,i,HZ_ALPHABETA,:),1),3),4))); hold on ; 
    h.FaceColor = CONTRAST_COLORS{i}; 
    errorbar(i*2,squeeze(mean(mean(mean(mdb(:,i,HZ_ALPHABETA,:),1),3),4)),squeeze(std(mean(mean(mdb(:,i,HZ_ALPHABETA,:),4),3),0,1))/sqrt(13),'k.'); 
    [h,p,ci,stats] = ttest(squeeze(mean(mean(mdb(:,1,HZ_ALPHABETA,:),3),4)),squeeze(mean(mean(mdb(:,i,HZ_ALPHABETA,:),3),4))); 
    if p < 0.05 ; text(i*2,-3.5,'*') ; end
end
set(gca,'XTick',2:2:10,'XTickLabel',CONTRAST_VALUES); xlabel('contrast(%)'); ylabel('dB'); title('alpha/beta, *p<0.05'); ylim([-3.75,0]);xlim([0,13]);

% peak frequency:
subplot(6,8,9); 
for i=1:5
    h = bar(i*2,squeeze(mean(peak_midgamma(:,i,1),1))); hold on ; 
    h.FaceColor = CONTRAST_COLORS{i}; 
    errorbar(i*2,squeeze(mean(peak_midgamma(:,i,1),1)),squeeze(std(peak_midgamma(:,i,1),0,1))/sqrt(13),'k.'); 
    [h,p,ci,stats] = ttest(peak_midgamma(:,1,1),peak_midgamma(:,i,1)); 
    if p < 0.05 ; text(i*2,55,'*') ; end
end
set(gca,'XTick',2:2:10,'XTickLabel',CONTRAST_VALUES); xlabel('contrast(%)'); ylabel('peak frequency(Hz)'); title('mid gamma, *p<0.05'); ylim([15,60]);xlim([0,13]);
subplot(6,8,10); 
for i=1:5
    h = bar(i*2,squeeze(mean(peak_lowgamma(:,i,1),1))); hold on ; 
    h.FaceColor = CONTRAST_COLORS{i}; 
    errorbar(i*2,squeeze(mean(peak_lowgamma(:,i,1),1)),squeeze(std(peak_lowgamma(:,i,1),0,1))/sqrt(13),'k.'); 
    [h,p,ci,stats] = ttest(peak_lowgamma(:,1,1),peak_lowgamma(:,i,1)); 
    if p < 0.05 ; text(i*2,34,'*') ; end
end
set(gca,'XTick',2:2:10,'XTickLabel',CONTRAST_VALUES); xlabel('contrast(%)'); ylabel('peak frequency(Hz)'); title('low gamma, *p<0.05'); ylim([15,35]);xlim([0,13]);
subplot(6,8,11); 
for i=1:5
    h = bar(i*2,squeeze(mean(peak_alpha(:,i,1),1))); hold on ; 
    h.FaceColor = CONTRAST_COLORS{i}; 
    errorbar(i*2,squeeze(mean(peak_alpha(:,i,1),1)),squeeze(std(peak_alpha(:,i,1),0,1))/sqrt(13),'k.'); 
    [h,p,ci,stats] = ttest(peak_alpha(:,1,1),peak_alpha(:,i,1)); 
    if p < 0.05 ; text(i*2,13.5,'*') ; end
end
set(gca,'XTick',2:2:10,'XTickLabel',CONTRAST_VALUES); xlabel('contrast(%)'); ylabel('peak frequency(Hz)'); title('alpha/beta, *p<0.05'); ylim([12,14]);xlim([0,13]);

subplot(2,3,4); 
for i=1:5
   shadedErrorBar(freqs(1:35),squeeze(mean(mean(mdb(:,i,1:35,:),1),4)),squeeze(std(mean(mdb(:,i,1:35,:),4),0,1))/sqrt(13),{'Color',CONTRAST_COLORS{i}}); hold on ;
   xlabel('frequency(hz)'); ylabel('dB'); xlim([0,70]); 
end
title('mean power across all orientations'); hline(0,'k'); vline([12.5,27.5,50]); ylim([-4,2]);
subplot(2,3,6) ; for i=1:5 ; plot(1,'Color',CONTRAST_COLORS{i},'LineWidth',5) ; hold on ; end ; legend({'100% contrast','50% contrast','25% contrast','15% contrast','5% contrast'});

clear ps; 
for i=1:5
    for j=1:60
        for k=1:180
            [h,p,ci,stats] = ttest(mdb(:,i,j,k)); 
            ps(i,j,k) = p; 
        end
    end
end

bps = ps < 0.01; 
for i=1:5
    total_ps = sum(sum(bps(i,[HZ_BOTHGAMMA,HZ_ALPHABETA],:))); 
    gamma_ps = sum(sum(bps(i,HZ_BOTHGAMMA,:))); 
    alpha_ps = sum(sum(bps(i,HZ_ALPHABETA,:)));     
    ratios(i) = gamma_ps/total_ps;
    ratios(i,2) = alpha_ps/total_ps;  
    ffrat(i) = gamma_ps/alpha_ps; 
end

figure,
for i=1:5
    subplot(6,8,i);
     imagesc(1:180,freqs,squeeze(ps(i,:,:)),[0,0.05]) ; axis xy ; colormap hot; title([num2str(CONTRAST_VALUES(i)), '% contrast']);  xlabel('orientation(deg\circ)'); 
     set(gca,'XTick',[45,90,135,180]) ; grid on; ylabel('frequency(hz)'); 
end
subplot(6,8,8) ; imagesc([0,0.05]) ;colormap hot; h = colorbar ; title(h,'p-value'); 

subplot(3,5,12); 
for i=1:5
    b = bar(i*2,ffrat(i)); b.FaceColor = CONTRAST_COLORS{i}; hold on;
end
xlabel('contrast(%)') ;set(gca,'XTick',2:2:10,'XTickLabel',CONTRAST_VALUES); ylabel('gamma/alpha ratio'); xlim([0,12]) ; title('FF/FB (gamma/alpha signficant dB)'); 

subplot(3,5,11); 
for i=1:5
   b = bar(i*3-2,ratios(i,1)); b.FaceColor = [0.1,0.15,0.2]; hold on ;
   b = bar(i*3-1,ratios(i,2)); b.FaceColor = [.5,0.45,0.4]; 
end
set(gca,'XTick',2:3:15,'XTickLabel',CONTRAST_VALUES); xlabel('contrast(%)') ;ylabel('significant dB (proportion)'); title('alpha/gamma significant dB')
subplot(3,5,15) ; plot(1,'Color',[0.1,0.15,0.2],'LineWidth',3); hold on ; plot(2,'Color',[.5,0.45,0.4],'LineWidth',3); legend({'FF (gamma)','FB (alpha)'}); 
subplot(3,5,14) ; for i=1:5 ; plot(1,'Color',CONTRAST_COLORS{i},'LineWidth',5) ; hold on ; end ; legend({'100% contrast','50% contrast','25% contrast','15% contrast','5% contrast'});


%{
for i=1:5
   b = bar(i*3-2,ratios(i,1)); b.FaceColor = [0.5,0,0]; hold on; 
   b = bar(i*3-1,ratios(i,2)); b.FaceColor = [0,0,0.5];    
end
set(gca,'XTick',1.5:3:15,'XTickLabel',CONTRAST_VALUES); xlabel('contrast(%)'); ylabel('
%}

% orientation selectivity indices: need single trials? 
cols1 = [1:90,90:-1:1];
cols2 = [1:45,45:-1:1,1:45,45:-1:1]; 
cols3 = [90:-1:1,1:90];
colsrgb(1,:,1) = mat2gray(cols1); colsrgb(1,:,2) = mat2gray(cols3); colsrgb(1,:,3) = mat2gray(cols2); 
figure,imagesc(colsrgb); set(gca,'YTick',[],'XTick',[45,90,135,180]) ; xlabel('orientation(deg\circ)'); 
figure,
%{
subplot(4,3,1);
imagesc(squeeze(mean(mean(mdb(:,:,[HZ_MIDGAMMA,HZ_LOWGAMMA],:),1),3))); colormap jet; h = colorbar ; title(h,'dB'); 
set(gca,'YTick',1:5,'YTickLabel',CONTRAST_VALUES); ylabel('contrast(%)'); xlabel('orientation(deg\circ)'); set(gca,'XTick',[45,90,135,180]); title('gamma amplitude (25-60Hz)'); hline(0.5:1:5.5,'k');

subplot(4,3,2); 
imagesc(squeeze(mean(mean(mdb(:,:,[HZ_ALPHABETA],:),1),3))); colormap jet; h = colorbar ; title(h,'dB'); 
set(gca,'YTick',1:5,'YTickLabel',CONTRAST_VALUES); ylabel('contrast(%)'); xlabel('orientation(deg\circ)'); set(gca,'XTick',[45,90,135,180]);  title('alpha/beta amplitdue (10-16Hz)');  hline(0.5:1:5.5,'k');
%}

subplot(4,4,1); 
for i=1:5
    plot(squeeze(mean(mean(mdb(:,i,[HZ_MIDGAMMA,HZ_LOWGAMMA],:),1),3))','Color',CONTRAST_COLORS{i},'LineWidth',2); hold on; 
end
title('gamma(25-60Hz)'); xlabel('orientation(deg\circ)'); ylabel('dB'); ylim([-.25,1.4]); set(gca,'XTick',[45,90,135,180]);  xlim([0,181]);

subplot(4,4,2); 
for i=1:5
    plot(squeeze(mean(mean(mdb(:,i,[HZ_ALPHABETA],:),1),3))','Color',CONTRAST_COLORS{i},'LineWidth',2); hold on; 
end
title('alpha/beta (10-16Hz)'); xlabel('orientation(deg\circ)'); ylabel('dB'); ylim([-3.5,-2]); set(gca,'XTick',[45,90,135,180]); xlim([0,181]);

subplot(4,4,3); 
plot(zscore(squeeze(mean(mean(mean(mdb(:,:,[HZ_MIDGAMMA,HZ_LOWGAMMA],:),1),2),3))),'Color',[0.1,0.15,0.2],'LineWidth',3); hold on ; 
plot(zscore(squeeze(mean(mean(mean(mdb(:,:,[HZ_ALPHABETA],:),1),2),3))),'Color',[.5,0.45,0.4],'LineWidth',3); ylim([-4,3]); xlabel('orientation(deg\circ)'); ylabel('z-sscore'); 
set(gca,'XTick',[45,90,135,180]); xlim([0,181]); 
%for i=1:length(oblique) ; text(oblique(i)-1.5,2.7,'\bullet','Color',squeeze(colsrgb(1,45,:))); end; 
%for i=1:length(vertical) ; text(vertical(i)-1.5,2.7,'\bullet','Color',squeeze(colsrgb(1,90,:))); end; 
%for i=1:length(horizontal) ; text(horizontal(i)-1.5,2.7,'\bullet','Color',squeeze(colsrgb(1,1,:))); end; 
%for i=1:length(cardinal) ; text(cardinal(i)-1.5,2.4,'\bullet','Color',squeeze(colsrgb(1,90,:))/2 + squeeze(colsrgb(1,1,:))/2); end; 
title('mean across all contrasts'); 

subplot(8,8,23) ; plot(1,'Color',[0.1,0.15,0.2],'LineWidth',3);hold on ; plot(1,'Color',[.5,0.45,0.4],'LineWidth',3); legend({'gamma (25-60Hz)','alpha/beta (10-16Hz)'});
subplot(8,8,24); for i=1:length(CONTRAST_COLORS) ; plot(1,'Color',CONTRAST_COLORS{i},'LineWidth',3); hold on; end ; legend({'100% contrast','50% contrast','25% contrast','15% contrast','5% contrast'}); 
subplot(8,8,32) ; imagesc(colsrgb); xlabel('orientation(deg\circ)'); set(gca,'XTick',[45,90,135,180],'XTickLabel',[45,90,135,180],'YTick',[]); 
subplot(8,8,31) ; plot(1,'Color',colsrgb(1,45,:),'LineWidth',3); hold on ;plot(1,'Color',colsrgb(1,1,:)/2 + colsrgb(1,90,:)/2,'LineWidth',3);  plot(1,'Color',colsrgb(1,90,:),'LineWidth',2) ; plot(1,'Color',colsrgb(1,1,:),'LineWidth',3) ; 
legend({'oblique (45\circ+/-10\circ, 135\circ+/-10\circ)','cardinal (90\circ+/-10\circ, 180\circ+/-10\circ)','vertical (90\circ+/-10\circ)','horizontal (180\circ+/-10\circ)'});

smooth_mdb = mdb; 
for i=1:size(mdb,1)
    for j=1:size(mdb,2)
        for k=1:size(mdb,3)
            smooth_mdb(i,j,k,:) = imfilter(squeeze(mdb(i,j,k,:)),fspecial('gaussian',[9,1],5)); 
        end
    end
end

HZ_ARRAY = {[HZ_MIDGAMMA,HZ_LOWGAMMA],HZ_ALPHABETA};
SORT_ARRAY = {'descend','ascend'};
clear hz_sorted allsubsi; 
for h=1:length(HZ_ARRAY)
mdb_hzi = squeeze(mean(mean(mdb(:,:,HZ_ARRAY{h},:),1),3)); 
mdb_hzi = imfilter(mdb_hzi,fspecial('gaussian',[1,9],5)); 
[sv,si] = sort(mdb_hzi,2,SORT_ARRAY{h}); 

[subsv,subsi] = sort(squeeze(mean(smooth_mdb(:,:,HZ_ARRAY{h},:),3)),3,SORT_ARRAY{h}); 

allsubsi(h,:,:,:) = subsi; 
allsi(h,:,:) = si; 
for i=1:5
   sortrgb(i,:,:) = colsrgb(1,si(i,:),:);  
end

hz_sorted(h,:,:,:) = sortrgb; 

end

m_gammadb = squeeze(mean(mean(mdb(:,:,[HZ_LOWGAMMA,HZ_MIDGAMMA],:),1),3)); 
plotinds = [1,3,5,7,9]; ylims = {[0,1.3],[0,1.3],[0,.6],[-.05,.4],[-.3,0]};
for i=1:5
    h=subplot(12,2,plotinds(i));
    [sv,si] = sort(m_gammadb(i,:),'descend'); 
    pos = get(h,'Position'); 
    pos(4) = pos(4) + 0.008;    
    set( h, 'Position', pos ) ;
    for j=1:length(sv)
        b = bar(j,sv(j)); hold on; 
        b.FaceColor = colsrgb(1,si(j),:); 
        b.EdgeColor = colsrgb(1,si(j),:);
    end
    plot(sv,'k'); xlim([0,180]); ylim(ylims{i}); ylabel([num2str(CONTRAST_VALUES(i)),'% (dB)']); 
    if i~=5
       set(gca,'XTick',[]); 
    else
        set(gca,'XTick',0:10:180); xlabel('sorted index (descending)');
    end
end

m_alphadb = squeeze(mean(mean(mdb(:,:,[HZ_ALPHABETA],:),1),3)); 
plotinds = [1,3,5,7,9]; ylims = {[0,1.3],[0,1.3],[0,.6],[-.05,.4],[-.3,0]};
for i=1:5
    h=subplot(12,2,plotinds(i));
    [sv,si] = sort(m_alphadb(i,:),'ascend'); 
    for j=1:length(sv)
        b = bar(j,sv(j)); hold on; 
        b.FaceColor = colsrgb(1,si(j),:); 
        b.EdgeColor = colsrgb(1,si(j),:);
    end
    pos = get(h,'Position'); 
    pos(4) = pos(4) + 0.008;    
    set( h, 'Position', pos ) ;
    plot(sv,'k'); xlim([0,180]); %title([num2str(CONTRAST_VALUES(i)),'% contrast']); 
    ylabel([num2str(CONTRAST_VALUES(i)),'% (dB)']); 
    if i~=5
       set(gca,'XTick',[]); 
    else
        set(gca,'XTick',0:10:180); xlabel('sorted index (ascending)');
    end
end

subplot(3,3,1); 
symbols = {'o','s','d','x','*'};
for c=1:5
[sv,si] = sort(m_gammadb(c,:),'descend'); 
for i=1:180
   plot(si(i),sv(i),symbols{c},'Color', squeeze(colsrgb(1,si(i),:)),'LineWidth',1); hold on; 
end
xlim([0,180]); ylim([-.3,1.4]); %plot(si,sv,'.','Color',CONTRAST_COLORS{c}); 
set(gca,'XTick',[45,90,135,180]) ; title('gamma (all orientation, all contrast)'); xlabel('orientation(deg\circ)'); ylabel('dB'); 
end

subplot(3,3,2); 
for c=1:5
[sv,si] = sort(m_alphadb(c,:),'descend'); 
for i=1:180
   plot(si(i),sv(i),symbols{c},'Color', squeeze(colsrgb(1,si(i),:)),'LineWidth',1); hold on; 
end
xlim([0,180]); ylim([-4,-2]); %plot(si,sv,'.','Color',CONTRAST_COLORS{c}); 
set(gca,'XTick',[45,90,135,180]) ; title('alpha (all orientation, all contrast)');  xlabel('orientation(deg\circ)'); ylabel('dB'); 
end

subplot(3,3,4); plot(1,'ko'); hold on; plot(1,'ks'); plot(1,'kd'); plot(1,'kx'); plot(1,'k*'); legend({'100%','50%','25%','15%','5%'});
subplot(18,6,60); imagesc(colsrgb); set(gca,'YTick',[],'XTick',[45,90,135,180]) ; xlabel('orientation(deg\circ)'); 

subplot(4,3,4) ;
imagesc(squeeze(hz_sorted(1,:,:,:))); 
xlabel('sorted index (gamma amplitude, \leftarrowhigher  lower\rightarrow)'); title('sorted gamma amplitude (25-60Hz)'); 
set(gca,'YTick',1:5,'YTickLabel',CONTRAST_VALUES);  ylabel('contrast(%)'); hline(0.5:1:5.5,'k'); 

subplot(4,3,5) ;
imagesc(squeeze(hz_sorted(2,:,:,:))); 
xlabel('sorted index (alpha/beta amplitude, \leftarrowlower  higher\rightarrow)'); title('sorted alpha/beta amplitude (10-16Hz)'); 
set(gca,'YTick',1:5,'YTickLabel',CONTRAST_VALUES); ylabel('contrast(%)'); hline(0.5:1:5.5,'k'); 

for i=1:5
    prop_oblique_gamma(i) = length(intersect(squeeze(allsi(1,i,1:60)),oblique))/60; 
    prop_cardinal_gamma(i) = length(intersect(squeeze(allsi(1,i,1:60)),cardinal))/60; 
    
    for j=1:13
        propsb_oblique_gamma(i,j) = length(intersect(squeeze(allsubsi(1,j,i,1:60)),oblique))/60; 
        propsb_cardinal_gamma(i,j) = length(intersect(squeeze(allsubsi(1,j,i,1:60)),cardinal))/60; 
        
        propsb_vertical_gamma(i,j) = length(intersect(squeeze(allsubsi(1,j,i,1:60)),vertical))/60; 
        propsb_horizontal_gamma(i,j) = length(intersect(squeeze(allsubsi(1,j,i,1:60)),horizontal))/60; 
        
        propsb_vertical_alphabeta(i,j) = length(intersect(squeeze(allsubsi(2,j,i,1:60)),vertical))/60; 
        propsb_horizontal_alphabeta(i,j) = length(intersect(squeeze(allsubsi(2,j,i,1:60)),horizontal))/60; 
        
        propsb_oblique_alphabeta(i,j) = length(intersect(squeeze(allsubsi(2,j,i,1:60)),oblique))/60; 
        propsb_cardinal_alphabeta(i,j) = length(intersect(squeeze(allsubsi(2,j,i,1:60)),cardinal))/60;         
    end      
end

subplot(4,3,7); 
for i=1:5
   h1 = bar(i*3-1,mean(propsb_oblique_gamma(i,:))); hold on; 
   h1.FaceColor = colsrgb(1,45,:);
   errorbar(i*3-1,mean(propsb_oblique_gamma(i,:)),std(propsb_oblique_gamma(i,:),0,2)/sqrt(13),'k.'); 
   h2 = bar(i*3,mean(propsb_cardinal_gamma(i,:))); 
   h2.FaceColor = (colsrgb(1,1,:) + colsrgb(1,90,:))/2; 
   errorbar(i*3,mean(propsb_cardinal_gamma(i,:)),std(propsb_cardinal_gamma(i,:),0,2)/sqrt(13),'k.'); 
   [h,ps(i),ci,stats] = ttest(propsb_oblique_gamma(i,:),propsb_cardinal_gamma(i,:)); 
   
   pmax = max([mean(propsb_oblique_gamma(i,:))+std(propsb_oblique_gamma(i,:),0,2)/sqrt(13),mean(propsb_cardinal_gamma(i,:))+std(propsb_cardinal_gamma(i,:),0,2)/sqrt(13)]); 
   text(i*3-1,pmax+pmax/10,format_p(ps(i))); 
   
end
ylim([0,.6]);
title('gamma oblique vs cardinal');  set(gca,'XTick',2.5:3:15,'XTickLabel',CONTRAST_VALUES);xlabel('contrast(%)');ylabel('proportion of top 1/3');

subplot(4,3,10); 
for i=1:5
   h1 = bar(i*3-1,mean(propsb_vertical_alphabeta(i,:))); hold on; 
   h1.FaceColor = colsrgb(1,90,:);
   errorbar(i*3-1,mean(propsb_vertical_alphabeta(i,:)),std(propsb_vertical_alphabeta(i,:),0,2)/sqrt(13),'k.');    
   h2 = bar(i*3,mean(propsb_horizontal_alphabeta(i,:))); 
   h2.FaceColor = (colsrgb(1,1,:) + colsrgb(1,180,:))/2; 
   errorbar(i*3,mean(propsb_horizontal_alphabeta(i,:)),std(propsb_horizontal_alphabeta(i,:),0,2)/sqrt(13),'k.');  
   
   [h,ps(i),ci,stats] = ttest(propsb_vertical_alphabeta(i,:),propsb_horizontal_alphabeta(i,:)); 
   
   pmax = max([mean(propsb_vertical_alphabeta(i,:))+std(propsb_vertical_alphabeta(i,:),0,2)/sqrt(13),mean(propsb_horizontal_alphabeta(i,:))+std(propsb_horizontal_alphabeta(i,:),0,2)/sqrt(13)]); 
   text(i*3-1,pmax+pmax/10,format_p(ps(i))); 
   
end
ylim([0,.35]); title('alpha/beta vertical vs horizontal'); set(gca,'XTick',2.5:3:15,'XTickLabel',CONTRAST_VALUES);xlabel('contrast(%)');ylabel('proportion of top 1/3');
subplot(4,6,17) ; plot(1,'Color',colsrgb(1,45,:),'LineWidth',3); hold on ;plot(1,'Color',(colsrgb(1,90,:) + colsrgb(1,1,:))/2,'LineWidth',3); legend({'oblique','cardinal'});
subplot(4,6,18) ; plot(1,'Color',colsrgb(1,90,:),'LineWidth',3); hold on ;plot(1,'Color',(colsrgb(1,1,:) + colsrgb(1,180,:))/2,'LineWidth',3); legend({'vertical','horizontal'});

subplot(4,3,8); 
for i=1:5
   h1 = bar(i*3-1,mean(propsb_vertical_gamma(i,:))); hold on; 
   h1.FaceColor = colsrgb(1,90,:);
   errorbar(i*3-1,mean(propsb_vertical_gamma(i,:)),std(propsb_vertical_gamma(i,:),0,2)/sqrt(13),'k.'); 
   h2 = bar(i*3,mean(propsb_horizontal_gamma(i,:))); 
   h2.FaceColor = (colsrgb(1,1,:) + colsrgb(1,180,:))/2; 
   errorbar(i*3,mean(propsb_horizontal_gamma(i,:)),std(propsb_horizontal_gamma(i,:),0,2)/sqrt(13),'k.'); 
   [h,ps(i),ci,stats] = ttest(propsb_vertical_gamma(i,:),propsb_horizontal_gamma(i,:)); 
   
   pmax = max([mean(propsb_vertical_gamma(i,:))+std(propsb_vertical_gamma(i,:),0,2)/sqrt(13),mean(propsb_horizontal_gamma(i,:))+std(propsb_horizontal_gamma(i,:),0,2)/sqrt(13)]); 
   text(i*3-1,pmax+pmax/10,format_p(ps(i))); 
   
end
ylim([0,.2]); ylabel('proportion of top 1/3');
title('gamma vertical vs horizontal');  set(gca,'XTick',2.5:3:15,'XTickLabel',CONTRAST_VALUES);xlabel('contrast(%)');

subplot(4,3,11); 
for i=1:5
   h1 = bar(i*3-1,mean(propsb_oblique_alphabeta(i,:))); hold on; 
   h1.FaceColor = colsrgb(1,45,:);
   errorbar(i*3-1,mean(propsb_oblique_alphabeta(i,:)),std(propsb_oblique_alphabeta(i,:),0,2)/sqrt(13),'k.');    
   h2 = bar(i*3,mean(propsb_cardinal_alphabeta(i,:))); 
   h2.FaceColor = (colsrgb(1,1,:) + colsrgb(1,90,:))/2; 
   errorbar(i*3,mean(propsb_cardinal_alphabeta(i,:)),std(propsb_cardinal_alphabeta(i,:),0,2)/sqrt(13),'k.');  
   
   [h,ps(i),ci,stats] = ttest(propsb_oblique_alphabeta(i,:),propsb_cardinal_alphabeta(i,:)); 
   
   pmax = max([mean(propsb_oblique_alphabeta(i,:))+std(propsb_oblique_alphabeta(i,:),0,2)/sqrt(13),mean(propsb_cardinal_alphabeta(i,:))+std(propsb_cardinal_alphabeta(i,:),0,2)/sqrt(13)]); 
   text(i*3-1,pmax+pmax/10,format_p(ps(i))); 
   
end
ylim([0,.45]); title('alpha/beta oblique vs cardinal'); set(gca,'XTick',2.5:3:15,'XTickLabel',CONTRAST_VALUES); xlabel('contrast(%)'); ylabel('proportion of top 1/3');


subplot(4,6,23) ; plot(1,'Color',colsrgb(1,45,:),'LineWidth',3); hold on ;plot(1,'Color',(colsrgb(1,90,:) + colsrgb(1,1,:))/2,'LineWidth',3); legend({'oblique','cardinal'});
subplot(4,6,24) ; plot(1,'Color',colsrgb(1,90,:),'LineWidth',3); hold on ;plot(1,'Color',(colsrgb(1,1,:) + colsrgb(1,180,:))/2,'LineWidth',3); legend({'vertical','horizontal'});


% FMRI RESULTS
cd e:/saved; 
allhigh = load('allhigh'); allhigh = allhigh.allhigh; 
alllow = load('alllow'); alllow = alllow.alllow; 

boldhigh = squeeze(mean(allhigh(2:end,:,:),1)); 
boldlow = squeeze(mean(alllow(2:end,:,:),1)); 
bold_contrast_labels = {'100','5'};
% figure 1a contrast orientation tuning curves
%figure 1a
figure,
subplot(3,3,1);
shadedErrorBar([],squeeze(mean(mean(allhigh(1:end,:,:),1),3)),squeeze(std(mean(allhigh(1:end,:,:),1),0,3))/sqrt(14),{'Color',CONTRAST_COLORS{1}}); hold on; 
shadedErrorBar([],squeeze(mean(mean(alllow(1:end,:,:),1),3)),squeeze(std(mean(alllow(1:end,:,:),1),0,3))/sqrt(14),{'Color',CONTRAST_COLORS{5}});
xlabel('orientation (deg\circ)'); ylabel('BOLD amplitude (a.u)'); grid on; set(gca,'XTick',[45,90,135,180]); title('BOLD orientation tuning');  xlim([0,180]); 

subplot(3,3,2); 
mhigh = squeeze(mean(mean(allhigh,1),3)); 
mlow = squeeze(mean(mean(alllow,1),3)); 
[sv,si] = sort(mhigh,'descend'); 
for i=1:180
    plot(i,mhigh(i),'o','Color',squeeze(colsrgb(1,i,:))); hold on; 
    plot(i,mlow(i),'*','Color',squeeze(colsrgb(1,i,:))); 
    grid on; 
end
xlim([0,180]) ; set(gca,'XTick',[45,90,135,180]); xlabel('orientation(deg\circ)'); 
subplot(4,3,3) ; plot(1,'Color',CONTRAST_COLORS{1},'LineWidth',3) ; hold on ; plot(1,'Color',CONTRAST_COLORS{5},'LineWidth',3); legend({'100% contrast','5% contrast'});
subplot(4,3,10); plot(1,'ko'); hold on ; plot(1,'k*') ; legend('100%','5%'); 


[sv_bh,si_bh] = sort(mean(boldhigh,2),'descend'); 
[sv_bl,si_bl] = sort(mean(boldlow,2),'descend'); 

%{
sortbold(1,:,:) = colsrgb(1,si_bh,:); 
sortbold(2,:,:) = colsrgb(1,si_bl,:); 

subplot(4,3,4) ; 
imagesc(sortbold) ; xlabel('sorted index (BOLD amplitude, \leftarrowhigher lower\rightarrow)'); ylabel('contrast(%)'); set(gca,'YTick',1:2,'YTickLabel',{100,5}); title('sorted BOLD amplitude'); hline(1.5,'k'); 
%}
h=subplot(6,2,1);
pos = get(h,'Position'); 
pos(4) = pos(4) + 0.008;    
set( h, 'Position', pos ) ;
for i=1:180   
    b = bar(i,sv_bh(i)); hold on; 
    b.FaceColor = colsrgb(1,si_bh(i),:); b.EdgeColor = colsrgb(1,si_bh(i),:); 
end
plot(sv_bh,'k'); 
ylim([-2.1,1.8]); xlim([0,180]); set(gca,'XTick',0:10:180,'XTickLabel',[]); ylabel('100% (BOLD A.U)'); 

h=subplot(6,2,3);
pos = get(h,'Position'); 
pos(4) = pos(4) + 0.008;    
set( h, 'Position', pos ) ;
for i=1:180   
    b = bar(i,sv_bl(i)); hold on; 
    b.FaceColor = colsrgb(1,si_bl(i),:); b.EdgeColor = colsrgb(1,si_bl(i),:);
end
plot(sv_bl,'k'); 
ylim([-2.1,1.8]); xlim([0,180]);set(gca,'Xtick',0:10:180); ylabel('5% (BOLD A.U)'); xlabel('sorted index (descending)'); 


for i=1:14
    [sv_bh,si_bh] = sort((boldhigh(:,i)),'descend'); 
    [sv_bl,si_bl] = sort((boldlow(:,i)),'descend'); 
    b_oblique_high(i) = length(intersect(oblique,si_bh(1:60)))/60; 
    b_oblique_low(i) = length(intersect(oblique,si_bl(1:60)))/60; 
    b_cardinal_high(i) = length(intersect(cardinal,si_bh(1:60)))/60; 
    b_cardinal_low(i) = length(intersect(cardinal,si_bl(1:60)))/60; 
    b_vertical_high(i) = length(intersect(vertical,si_bh(1:60)))/60; 
    b_vertical_low(i) = length(intersect(vertical,si_bl(1:60)))/60; 
    b_horizontal_high(i) = length(intersect(horizontal,si_bh(1:60)))/60; 
    b_horizontal_low(i) = length(intersect(horizontal,si_bl(1:60)))/60; 
end

subplot(4,3,4); 
b = bar(1,mean(b_oblique_high)); hold on ; errorbar(1,mean(b_oblique_high),std(b_oblique_high)/sqrt(13),'k.'); 
b.FaceColor = colsrgb(1,45,:); 
b = bar(2,mean(b_cardinal_high)); errorbar(2,mean(b_cardinal_high),std(b_cardinal_high)/sqrt(13),'k.'); 
b.FaceColor = colsrgb(1,90,:)/2 + colsrgb(1,180,:)/2; 
ymax = max([mean(b_oblique_high)+std(b_oblique_high)/sqrt(13),mean(b_cardinal_high),std(b_cardinal_high)/sqrt(13)]);
[h,p,ci,stats] = ttest(b_oblique_high,b_cardinal_high); text(1,ymax+ymax/10,format_p(p)); 

b = bar(4,mean(b_oblique_low)); hold on ; errorbar(4,mean(b_oblique_low),std(b_oblique_low)/sqrt(13),'k.'); 
b.FaceColor = colsrgb(1,45,:); 
b = bar(5,mean(b_cardinal_low)); errorbar(5,mean(b_cardinal_low),std(b_cardinal_low)/sqrt(13),'k.'); 
b.FaceColor = colsrgb(1,90,:)/2 + colsrgb(1,180,:)/2; 
ymax = max([mean(b_oblique_low)+std(b_oblique_low)/sqrt(13),mean(b_cardinal_low),std(b_cardinal_low)/sqrt(13)]);
[h,p,ci,stats] = ttest(b_oblique_low,b_cardinal_low); text(4,ymax+ymax/10,format_p(p)); 

ylabel('proportion of top 1/3');set(gca,'XTick',[1.5,4.5],'XTickLabel',{100,5}); xlabel('contrast(%)'); title('BOLD oblique vs cardinal'); 

subplot(4,3,5); 
b = bar(1,mean(b_vertical_high)); hold on ; errorbar(1,mean(b_vertical_high),std(b_vertical_high)/sqrt(13),'k.'); 
b.FaceColor = colsrgb(1,90,:); 
b = bar(2,mean(b_horizontal_high)); errorbar(2,mean(b_horizontal_high),std(b_horizontal_high)/sqrt(13),'k.'); 
b.FaceColor = colsrgb(1,1,:);
ymax = max([mean(b_vertical_high)+std(b_vertical_high)/sqrt(13),mean(b_horizontal_high),std(b_horizontal_high)/sqrt(13)]);
[h,p,ci,stats] = ttest(b_vertical_high,b_horizontal_high); text(1,ymax+ymax/10,format_p(p)); 

b = bar(4,mean(b_vertical_low)); hold on ; errorbar(4,mean(b_vertical_low),std(b_vertical_low)/sqrt(13),'k.'); 
b.FaceColor = colsrgb(1,90,:); 
b = bar(5,mean(b_horizontal_low)); errorbar(5,mean(b_horizontal_low),std(b_horizontal_low)/sqrt(13),'k.'); 
b.FaceColor = colsrgb(1,1,:);
ylabel('proportion of top 1/3'); set(gca,'XTick',[1.5,4.5],'XTickLabel',{100,5}); xlabel('contrast(%)'); title('BOLD vertical vs horizontal'); 
ymax = max([mean(b_vertical_low)+std(b_vertical_low)/sqrt(13),mean(b_horizontal_low),std(b_horizontal_low)/sqrt(13)]);
[h,p,ci,stats] = ttest(b_vertical_low,b_horizontal_low); text(4,ymax+ymax/14,format_p(p)); 

subplot(8,4,25) ; plot(1,'Color',colsrgb(1,45,:),'LineWidth',3); hold on ;plot(1,'Color',(colsrgb(1,90,:) + colsrgb(1,1,:))/2,'LineWidth',3); legend({'oblique','cardinal'});
subplot(8,4,26) ; plot(1,'Color',colsrgb(1,90,:),'LineWidth',3); hold on ;plot(1,'Color',(colsrgb(1,1,:) + colsrgb(1,180,:))/2,'LineWidth',3); legend({'vertical','horizontal'});


subplot(8,4,31) ; imagesc(colsrgb); xlabel('orientation(deg\circ)'); set(gca,'XTick',[45,90,135,180],'XTickLabel',[45,90,135,180],'YTick',[]); 
subplot(8,4,32) ; plot(1,'Color',colsrgb(1,45,:),'LineWidth',3); hold on ;plot(1,'Color',colsrgb(1,1,:)/2 + colsrgb(1,90,:)/2,'LineWidth',3);  plot(1,'Color',colsrgb(1,90,:),'LineWidth',2) ; plot(1,'Color',colsrgb(1,1,:),'LineWidth',3) ; 
legend({'oblique (45\circ+/-10\circ, 135\circ+/-10\circ)','cardinal (90\circ+/-10\circ, 180\circ+/-10\circ)','vertical (90\circ+/-10\circ)','horizontal (180\circ+/-10\circ)'});


% COMBINING BOTH MODALITIES: FIGURE 4 

figure
clear highcorrs lowcorrs highps lowps
for i=1:5 
    for j=1:60
        [highcorrs(i,j),highps(i,j)] = corr(mean(boldhigh,2),squeeze(mean(mdb(:,i,j,:),1))); 
        [lowcorrs(i,j),lowps(i,j)] = corr(mean(boldlow,2),squeeze(mean(mdb(:,i,j,:),1))); 
    end
end

oblique_color = squeeze(colsrgb(1,45,:));
vertical_color = squeeze(colsrgb(1,90,:)); 
horizontal_color = squeeze(colsrgb(1,1,:)); 

subplot(4,3,1); 
plot(zscore(squeeze(mean(mean(allhigh(2:end,:,:),1),3)))); hold on ; 
plot(zscore(squeeze(mean(mean(mdb(:,1,[HZ_MIDGAMMA,HZ_LOWGAMMA],:),1),3)))); 
legend({'BOLD','EEG gamma'});xlabel('orientation(deg\circ)'); ylabel('tuning z-score'); xlim([0,180]) ; set(gca,'XTick',[45,90,135]); 
title('100% contrast'); 
subplot(4,3,2); 
plot(zscore(squeeze(mean(mean(alllow(2:end,:,:),1),3)))); hold on ; 
plot(zscore(squeeze(mean(mean(mdb(:,5,HZ_ALPHABETA,:),1),3)))); 
legend({'BOLD','EEG alpha/beta'});xlabel('orientation(deg\circ)'); ylabel('tuning z-score'); xlim([0,180]) ; set(gca,'XTick',[45,90,135]); 
title('5% contrast'); 

orientations = zeros(1,180) ; orientations([oblique,vertical,horizontal]) = 1; other_orientations = find(orientations==0); 

subplot(4,6,7); 
plot(squeeze(mean(mean(mdb(:,1,[HZ_MIDGAMMA,HZ_LOWGAMMA],other_orientations),1),3)),squeeze(mean(boldhigh(other_orientations,:),2)),'d','LineWidth',1,'Color',[0.3,0.3,0.3]);  hold on ; 
plot(squeeze(mean(mean(mdb(:,1,[HZ_MIDGAMMA,HZ_LOWGAMMA],oblique),1),3)),squeeze(mean(boldhigh(oblique,:),2)),'d','LineWidth',2,'Color',oblique_color);
plot(squeeze(mean(mean(mdb(:,1,[HZ_MIDGAMMA,HZ_LOWGAMMA],vertical),1),3)),squeeze(mean(boldhigh(vertical,:),2)),'d','LineWidth',2,'Color',vertical_color); 
plot(squeeze(mean(mean(mdb(:,1,[HZ_MIDGAMMA,HZ_LOWGAMMA],horizontal),1),3)),squeeze(mean(boldhigh(horizontal,:),2)),'d','LineWidth',2,'Color',horizontal_color); 
x = double(squeeze(mean(mean(mdb(:,1,[HZ_MIDGAMMA,HZ_LOWGAMMA],:),1),3))); y = double(squeeze(mean(boldhigh,2))); 
B = [ones(size(x(:)))  x(:)]\y(:);
y_fit = [ones(size(x(:)))  x(:)]*B;
plot(x,y_fit,'k','LineWidth',3); xlabel('gamma 100% contrast (dB)'); ylabel('BOLD 100% contrast (AU)'); title('100%contrast (gamma vs BOLD)'); 
[rho,p] = corr(x,y); text(min(x)-abs(min(x))/10,max(y)+abs(max(y))/10,format_rho(rho)); xlim([min(x)-abs(min(x))/10,max(x)+abs(max(x))/10]); ylim([min(y)-abs(min(y))/10,max(y)+abs(max(y))/4]); 

subplot(4,6,8); 
plot(squeeze(mean(mean(mdb(:,1,HZ_ALPHABETA,other_orientations),1),3)),squeeze(mean(boldhigh(other_orientations,:),2)),'d','LineWidth',1,'Color',[0.3,0.3,0.3]);  hold on ; 
plot(squeeze(mean(mean(mdb(:,1,HZ_ALPHABETA,oblique),1),3)),squeeze(mean(boldhigh(oblique,:),2)),'d','LineWidth',2,'Color',oblique_color);
plot(squeeze(mean(mean(mdb(:,1,HZ_ALPHABETA,vertical),1),3)),squeeze(mean(boldhigh(vertical,:),2)),'d','LineWidth',2,'Color',vertical_color); 
plot(squeeze(mean(mean(mdb(:,1,HZ_ALPHABETA,horizontal),1),3)),squeeze(mean(boldhigh(horizontal,:),2)),'d','LineWidth',2,'Color',horizontal_color); 
x = double(squeeze(mean(mean(mdb(:,1,HZ_ALPHABETA,:),1),3))); y = double(squeeze(mean(boldhigh,2))); 
B = [ones(size(x(:)))  x(:)]\y(:);
y_fit = [ones(size(x(:)))  x(:)]*B;
plot(x,y_fit,'k','LineWidth',3); xlabel('alpha/beta 100% contrast (dB)'); ylabel('BOLD 100% contrast (AU)'); title('100%contrast (alpha/beta vs BOLD)'); 
[rho,p] = corr(x,y); text(min(x)-abs(min(x))/10,max(y)+abs(max(y))/10,format_rho(rho));   xlim([min(x)-abs(min(x))/10,max(x)+abs(max(x))/10]); ylim([min(y)-abs(min(y))/10,max(y)+abs(max(y))/4]); 

subplot(4,6,9); 
plot(squeeze(mean(mean(mdb(:,5,[HZ_MIDGAMMA,HZ_LOWGAMMA],other_orientations),1),3)),squeeze(mean(boldlow(other_orientations,:),2)),'d','LineWidth',1,'Color',[0.3,0.3,0.3]);  hold on ; 
plot(squeeze(mean(mean(mdb(:,5,[HZ_MIDGAMMA,HZ_LOWGAMMA],oblique),1),3)),squeeze(mean(boldlow(oblique,:),2)),'d','LineWidth',2,'Color',oblique_color);
plot(squeeze(mean(mean(mdb(:,5,[HZ_MIDGAMMA,HZ_LOWGAMMA],vertical),1),3)),squeeze(mean(boldlow(vertical,:),2)),'d','LineWidth',2,'Color',vertical_color); 
plot(squeeze(mean(mean(mdb(:,5,[HZ_MIDGAMMA,HZ_LOWGAMMA],horizontal),1),3)),squeeze(mean(boldlow(horizontal,:),2)),'d','LineWidth',2,'Color',horizontal_color); 
x = double(squeeze(mean(mean(mdb(:,5,[HZ_MIDGAMMA,HZ_LOWGAMMA],:),1),3))); y = double(squeeze(mean(boldlow,2))); 
B = [ones(size(x(:)))  x(:)]\y(:);
y_fit = [ones(size(x(:)))  x(:)]*B;
plot(x,y_fit,'k','LineWidth',3); xlabel('gamma 5% contrast (dB)'); ylabel('BOLD 5% contrast (AU)'); title('5%contrast (gamma vs BOLD)'); 
[rho,p] = corr(x,y); text(min(x)-abs(min(x))/10,max(y)+abs(max(y))/10,format_rho(rho));  xlim([min(x)-abs(min(x))/10,max(x)+abs(max(x))/10]); ylim([min(y)-abs(min(y))/10,max(y)+abs(max(y))/4]); 

subplot(4,6,10); 
plot(squeeze(mean(mean(mdb(:,5,HZ_ALPHABETA,other_orientations),1),3)),squeeze(mean(boldlow(other_orientations,:),2)),'d','LineWidth',1,'Color',[0.3,0.3,0.3]);  hold on ; 
plot(squeeze(mean(mean(mdb(:,5,HZ_ALPHABETA,oblique),1),3)),squeeze(mean(boldlow(oblique,:),2)),'d','LineWidth',2,'Color',oblique_color);
plot(squeeze(mean(mean(mdb(:,5,HZ_ALPHABETA,vertical),1),3)),squeeze(mean(boldlow(vertical,:),2)),'d','LineWidth',2,'Color',vertical_color); 
plot(squeeze(mean(mean(mdb(:,5,HZ_ALPHABETA,horizontal),1),3)),squeeze(mean(boldlow(horizontal,:),2)),'d','LineWidth',2,'Color',horizontal_color); 
x = double(squeeze(mean(mean(mdb(:,5,HZ_ALPHABETA,:),1),3))); y = double(squeeze(mean(boldlow,2))); 
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
plot(freqs,squeeze(highcorrs(1,:)),'Color',CONTRAST_COLORS{1},'LineWidth',2); hold on ; hline(0,'k'); 
plot(freqs,squeeze(lowcorrs(5,:)),'Color',CONTRAST_COLORS{5},'LineWidth',2); hold on ; hline(0,'k'); 
xlim([0,120]); ylabel('rho'); xlabel('frequency(hz)'); title('EEG-BOLD tuning correlation'); 

subplot(4,3,8) ; 
plot(freqs,squeeze(highcorrs(1,:)).^2,'Color',CONTRAST_COLORS{1},'LineWidth',2); hold on ; hline(0,'k'); 
plot(freqs,squeeze(lowcorrs(5,:)).^2,'Color',CONTRAST_COLORS{5},'LineWidth',2); hold on ; hline(0,'k'); 
xlim([0,120]); ylabel('rho^{2}'); xlabel('frequency(hz)'); title('BOLD tuning variance explained'); 

subplot(4,3,9); plot(1,'Color',CONTRAST_COLORS{1},'LineWidth',5); hold on; plot(1,'Color',CONTRAST_COLORS{5},'LineWidth',5); legend({'100% contrast','5% contrast'});
subplot(4,3,12); plot(1,'Color',CONTRAST_COLORS{1},'LineWidth',5); hold on; plot(1,'Color',CONTRAST_COLORS{5},'LineWidth',5); legend({'100% contrast','5% contrast'});

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



%%%% figure for linearity of alpha/beta and gamma

figure,

gamma_100 = squeeze(mean(mean(mdb(:,1,HZ_MIDGAMMA,:),3),4)); 
gamma_50 = squeeze(mean(mean(mdb(:,2,HZ_MIDGAMMA,:),3),4))*2; 
gamma_25 = squeeze(mean(mean(mdb(:,3,HZ_MIDGAMMA,:),3),4))*4; 
gamma_15 = squeeze(mean(mean(mdb(:,4,HZ_MIDGAMMA,:),3),4))*6.6667; 
gamma_5 = squeeze(mean(mean(mdb(:,5,HZ_MIDGAMMA,:),3),4))*20; 

alphabeta_100 = squeeze(mean(mean(mdb(:,1,HZ_ALPHABETA,:),3),4)); 
alphabeta_50 = squeeze(mean(mean(mdb(:,2,HZ_ALPHABETA,:),3),4))*2; 
alphabeta_25 = squeeze(mean(mean(mdb(:,3,HZ_ALPHABETA,:),3),4))*4; 
alphabeta_15 = squeeze(mean(mean(mdb(:,4,HZ_ALPHABETA,:),3),4))*6.6667; 
alphabeta_5 = squeeze(mean(mean(mdb(:,5,HZ_ALPHABETA,:),3),4))*20; 


mspec_100 = squeeze(mean(mdb(:,1,:,:),4))*1; 
mspec_50 = squeeze(mean(mdb(:,2,:,:),4))*2; 
mspec_25 = squeeze(mean(mdb(:,3,:,:),4))*4; 
mspec_15 = squeeze(mean(mdb(:,4,:,:),4))*6.6667; 

subplot(4,4,2); 
shadedErrorBar(1:2:24,mean(mspec_100(:,1:12),1),std(mspec_100(:,1:12),0,1)/sqrt(13),{'Color',CONTRAST_COLORS{1}}); hold on; 
shadedErrorBar(1:2:24,mean(mspec_50(:,1:12),1),std(mspec_50(:,1:12),0,1)/sqrt(13),{'Color',CONTRAST_COLORS{2}}); 
shadedErrorBar(1:2:24,mean(mspec_25(:,1:12),1),std(mspec_25(:,1:12),0,1)/sqrt(13),{'Color',CONTRAST_COLORS{3}}); 
shadedErrorBar(1:2:24,mean(mspec_15(:,1:12),1),std(mspec_15(:,1:12),0,1)/sqrt(13),{'Color',CONTRAST_COLORS{4}}); xlim([0.5,23.5]); xlabel('frequency(hz)'); ylabel('dB'); 
title('alpha/beta'); 

subplot(4,4,1); 
shadedErrorBar(26:2:100,mean(mspec_100(:,13:50),1),std(mspec_100(:,13:50),0,1)/sqrt(13),{'Color',CONTRAST_COLORS{1}}); hold on; 
shadedErrorBar(26:2:100,mean(mspec_50(:,13:50),1),std(mspec_50(:,13:50),0,1)/sqrt(13),{'Color',CONTRAST_COLORS{2}}); 
shadedErrorBar(26:2:100,mean(mspec_25(:,13:50),1),std(mspec_25(:,13:50),0,1)/sqrt(13),{'Color',CONTRAST_COLORS{3}}); 
shadedErrorBar(26:2:100,mean(mspec_15(:,13:50),1),std(mspec_15(:,13:50),0,1)/sqrt(13),{'Color',CONTRAST_COLORS{4}}); xlim([25,101]); xlabel('frequency(hz)'); ylabel('dB'); 
title('gamma'); 
subplot(4,4,4); 
plot(1,'Color',CONTRAST_COLORS{1},'LineWidth',3); hold on; plot(1,'Color',CONTRAST_COLORS{2},'LineWidth',3);plot(1,'Color',CONTRAST_COLORS{3},'LineWidth',3);plot(1,'Color',CONTRAST_COLORS{4},'LineWidth',3);
legend({'100%','50% (x2)','25% (x4)','15% (x6.67)'});

subplot(4,4,5) ;
h = bar(1,mean(gamma_100)); hold on; errorbar(1,mean(gamma_100),std(gamma_100)/sqrt(13),'k.'); 
h.FaceColor = CONTRAST_COLORS{1}; 

h = bar(3,mean(gamma_50));errorbar(3,mean(gamma_50),std(gamma_50)/sqrt(13),'k.'); 
h.FaceColor = CONTRAST_COLORS{2}; 
[h,p,ci,stats] = ttest(gamma_100,gamma_50); maxval = mean(gamma_50) + std(gamma_50)/sqrt(13); 
text(3,maxval+maxval/10,format_p(p)); 

h = bar(5,mean(gamma_25)); errorbar(5,mean(gamma_25),std(gamma_25)/sqrt(13),'k.'); 
h.FaceColor = CONTRAST_COLORS{3}; 
[h,p,ci,stats] = ttest(gamma_100,gamma_25); maxval = mean(gamma_25) + std(gamma_25)/sqrt(13); 
text(5,maxval+maxval/10,format_p(p)); 

h = bar(7,mean(gamma_15)); errorbar(7,mean(gamma_15),std(gamma_15)/sqrt(13),'k.'); 
h.FaceColor = CONTRAST_COLORS{4}; 
[h,p,ci,stats] = ttest(gamma_100,gamma_15); maxval = mean(gamma_15) + std(gamma_15)/sqrt(13); 
text(7,maxval+maxval/10,format_p(p)); 

set(gca,'XTick',1:2:7,'XTickLabel',{'100','50 (x2)','25 (x4)','15 (x6.67)'});xlabel('contrast(%)'); ylabel('dB'); 
title('gamma (40-60Hz)');

subplot(4,4,6) ;
h = bar(1,mean(alphabeta_100)); hold on; errorbar(1,mean(alphabeta_100),std(alphabeta_100)/sqrt(13),'k.'); 
h.FaceColor = CONTRAST_COLORS{1};

h = bar(3,mean(alphabeta_50));errorbar(3,mean(alphabeta_50),std(alphabeta_50)/sqrt(13),'k.'); 
h.FaceColor = CONTRAST_COLORS{2}; 
[h,p,ci,stats] = ttest(alphabeta_100,alphabeta_50); 
text(3,2,format_p(p)); 

h = bar(5,mean(alphabeta_25)); errorbar(5,mean(alphabeta_25),std(alphabeta_25)/sqrt(13),'k.'); 
h.FaceColor = CONTRAST_COLORS{3}; 
[h,p,ci,stats] = ttest(alphabeta_100,alphabeta_25); 
text(5,2,format_p(p)); 

h = bar(7,mean(alphabeta_15)); errorbar(7,mean(alphabeta_15),std(alphabeta_15)/sqrt(13),'k.'); 
h.FaceColor = CONTRAST_COLORS{4}; 
[h,p,ci,stats] = ttest(alphabeta_100,alphabeta_15); 
text(7,2,format_p(p)); 

set(gca,'XTick',1:2:7,'XTickLabel',{'100','50 (x2)','25 (x4)','15 (x6.67)'});xlabel('contrast(%)'); ylabel('dB'); 
title('alpha/beta (10-16Hz)'); ylim([-25,5]);

cd C:\Users\butr2901\Desktop\tuning2

cd e:/saved; sub_roits = load('sub_roits'); sub_roits = sub_roits.sub_roits; 
for i=1:9
   for j=1:60
      highcorrs(i,j) = corr(squeeze(mean(mdb(:,1,j,:))),squeeze(mean(sub_roits(:,i,1,:),1)));  
      lowcorrs(i,j) = corr(squeeze(mean(mdb(:,5,j,:))),squeeze(mean(sub_roits(:,i,2,:),1)));  
   end   
end
roi_colors = {[153,32,102]/255,[0,0,255]/255,[0,153,255]/255,[0,255,255]/255,[0,255,0]/255,[255,255,0]/255,[255,105,0]/255,[255,0,0]/255,[204,16,51]/255}; 
rnames = {'V1','V2','V3','hV4','VO','PH','MST','hMT','LO'};

for i=1:9
   plot(lowcorrs(i,:),'Color',roi_colors{i},'LineWidth',2);  hold on; 
end

subplot(1,3,1) ; imagesc(freqs,1:9,highcorrs,[-1,1]); xlabel('frequency(hz)');set(gca,'YTick',1:9,'YTickLabel',rnames); colormap jet; title('100% contrast'); 
subplot(1,3,2); imagesc(freqs,1:9,lowcorrs,[-1,1]); xlabel('frequency(hz)'); set(gca,'YTick',1:9,'YTickLabel',rnames); colormap jet; title('5% contrast'); 
subplot(1,3,3) ; imagesc([-1,1]) ; h = colorbar; title(h,'rho'); 












