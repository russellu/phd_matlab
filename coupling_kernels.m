clear all  ;close all ; 
networks = {'dmn','anterior-vis','lateral-vis','motor','auditory','right-occ','left-occ','left-exec','right-exec','alien'} ;
cd C:\shared\savecoupling
comps = load_untouch_nii('comps.nii.gz') ; 
gamma1 = load('gamma_01') ; gamma1 = gamma1.allxcorrs ;
gamma2 = load('gamma_02') ; gamma2 = gamma2.allxcorrs ; 
movie = load('movie') ; movie = movie.allxcorrs ; 
rest = load('rest') ; rest = rest.allxcorrs ; 
mgamma = (gamma1 + gamma2) / 2 ; 
gamma = squeeze(mean(mgamma(:,:,:,1:200,:),4)) ; 
movie = squeeze(mean(movie(:,:,:,1:200,:),4)) ;
rest = squeeze(mean(rest(:,:,:,1:200,:),4)) ;
freqs = 1:2:90 ; times = (-20:20)*0.693 ; 
subplot(2,2,1) ; imagesc(times,freqs,squeeze(mean(mean(gamma(:,6:7,:,:),1),2)),[-.2,.2]) ; axis xy ; vline(0,'k') ; ylabel('frequency(hz)') ; xlabel('time(s)') ; 
subplot(2,2,2) ; imagesc(times,freqs,squeeze(mean(mean(rest(:,6:7,:,:),1),2)),[-.2,.2]) ; axis xy ; vline(0,'k') ; ylabel('frequency(hz)') ; xlabel('time(s)') ; 
subplot(2,2,3) ; imagesc(times,freqs,squeeze(mean(mean(movie(:,6:7,:,:),1),2)),[-.2,.2]) ; axis xy ; vline(0,'k') ; ylabel('frequency(hz)') ; xlabel('time(s)') ; 
figure,subplot(2,2,1) ; imagesc(times,freqs,squeeze(mean(mean(gamma(:,1,:,:),1),2)),[-.2,.2]) ; axis xy ; vline(0,'k') ; ylabel('frequency(hz)') ; xlabel('time(s)') ; 
subplot(2,2,2) ; imagesc(times,freqs,squeeze(mean(mean(rest(:,1,:,:),1),2)),[-.2,.2]) ; axis xy ; vline(0,'k') ; ylabel('frequency(hz)') ; xlabel('time(s)') ; 
subplot(2,2,3) ; imagesc(times,freqs,squeeze(mean(mean(movie(:,1,:,:),1),2)),[-.2,.2]) ; axis xy ; vline(0,'k') ; ylabel('frequency(hz)') ; xlabel('time(s)') ; 
figure,subplot(2,2,1) ; imagesc(times,freqs,squeeze(mean(mean(gamma(:,2,:,:),1),2)),[-.2,.2]) ; axis xy ; vline(0,'k') ; ylabel('frequency(hz)') ; xlabel('time(s)') ; 
subplot(2,2,2) ; imagesc(times,freqs,squeeze(mean(mean(rest(:,2,:,:),1),2)),[-.2,.2]) ; axis xy ; vline(0,'k') ; ylabel('frequency(hz)') ; xlabel('time(s)') ; 
subplot(2,2,3) ; imagesc(times,freqs,squeeze(mean(mean(movie(:,2,:,:),1),2)),[-.2,.2]) ; axis xy ; vline(0,'k') ; ylabel('frequency(hz)') ; xlabel('time(s)') ; 

hrf = spm_hrf(0.693) ; 
trf = zeros(1,41) ; trf(21:end) = hrf(1:21) ; 

subplot(2,2,1) ; 
shadedErrorBar([],squeeze(mean(mean(mean(rest(:,6:7,4:8,:),1),2),3)),squeeze(std(mean(mean(rest(:,6:7,4:8,:),2),3),0,1))./sqrt(7),{'b'}) ; hold on ; 
shadedErrorBar([],squeeze(mean(mean(mean(movie(:,6:7,4:8,:),1),2),3)),squeeze(std(mean(mean(movie(:,6:7,4:8,:),2),3),0,1))./sqrt(7),{'r'}) ; 
plot(trf,'k','LineWidth',3) ; ylim([-.28,.17]) ; 
xlim([1,41]) ; hline(0,'k') ; set(gca,'XTick',1:5:length(times),'XTickLabel',round(times(1:5:end)*10)/10) ; vline(21,'r') ; 
ylabel('correlation(rho)') ; xlabel('time lag (s)') ; 
subplot(2,2,2) ; 
shadedErrorBar([],squeeze(mean(mean(mean(rest(:,1,8:14,:),1),2),3)),squeeze(std(mean(mean(rest(:,1,8:14,:),2),3),0,1))./sqrt(7),{'b'}) ; hold on ; 
shadedErrorBar([],squeeze(mean(mean(mean(movie(:,1,8:14,:),1),2),3)),squeeze(std(mean(mean(movie(:,1,8:14,:),2),3),0,1))./sqrt(7),{'r'}) ; 
plot(trf,'k','LineWidth',3) ; ylim([-.28,.2]) ; 
xlim([1,41]) ; hline(0,'k') ; set(gca,'XTick',1:5:length(times),'XTickLabel',round(times(1:5:end)*10)/10) ; vline(21,'r') ; 
ylabel('correlation(rho)') ; xlabel('time lag (s)') ; 
subplot(2,2,3) ; 
shadedErrorBar([],squeeze(mean(mean(mean(gamma(:,6:7,20:30,:),1),2),3)),squeeze(std(mean(mean(gamma(:,6:7,20:30,:),2),3),0,1))./sqrt(7),{'m'}) ; hold on ; 
shadedErrorBar([],squeeze(mean(mean(mean(rest(:,6:7,20:30,:),1),2),3)),squeeze(std(mean(mean(rest(:,6:7,20:30,:),2),3),0,1))./sqrt(7),{'b'}) ; 
shadedErrorBar([],squeeze(mean(mean(mean(movie(:,6:7,20:30,:),1),2),3)),squeeze(std(mean(mean(movie(:,6:7,20:30,:),2),3),0,1))./sqrt(7),{'r'}) ; 
plot(trf,'k','LineWidth',3) ; ylim([-.12,.16]) ; 
xlim([1,41]) ; hline(0,'k') ; set(gca,'XTick',1:5:length(times),'XTickLabel',round(times(1:5:end)*10)/10) ; vline(21,'r') ; 
ylabel('correlation(rho)') ; xlabel('time lag (s)') ; 



%{
anat = load_untouch_nii('mean.nii') ; 
dmn = comps.img(:,:,:,1) ; occ = comps.img(:,:,:,6) + comps.img(:,:,:,7) ; ling = comps.img(:,:,:,2) ; 
subplot(1,3,1) ; 
plotoverlayIntensity2D(anat.img(:,:,15),mat2gray(mean(dmn,3)),mean(dmn,3),270) ; 
subplot(1,3,2) ; 
plotoverlayIntensity2D(anat.img(:,:,10),mat2gray(mean(occ,3)),mean(occ,3),270) ; 
subplot(1,3,3) ; 
plotoverlayIntensity2D(anat.img(:,:,12),mat2gray(mean(ling,3)),mean(ling,3),270) ; 
%}











