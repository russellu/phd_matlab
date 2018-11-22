clear all ; close all ;
rsns = {'anterior_visual.nii','dmn.nii','higher_visual.nii','primary_visual.nii','motor.nii'};
rsn_titles = {'Anterior Visual','Default Mode','Higher-order Visual','Primary Visual','Motor'};
cd E:\sim\saved
% (sub,scan,freq,rsn,time)
times = (-25:25)*0.693; 
yes_gsr = load('mc_yes_gsr_xcorrs.mat'); yes_gsr = yes_gsr.mc_yes_gsr_xcorrs; 
no_gsr = load('mc_no_gsr_xcorrs.mat'); no_gsr = no_gsr.mc_no_gsr_xcorrs; 

temp = zeros(size(yes_gsr,1),4,size(yes_gsr,3),size(yes_gsr,4),size(yes_gsr,5));
temp(:,1,:,:,:) = squeeze(mean(yes_gsr(:,1:2,:,:,:),2)); 
temp(:,2,:,:,:) = squeeze(mean(yes_gsr(:,3:4,:,:,:),2)); 
temp(:,3,:,:,:) = squeeze(mean(yes_gsr(:,5,:,:,:),2)); 
temp(:,4,:,:,:) = squeeze(mean(yes_gsr(:,6,:,:,:),2)); 
yes_gsr = temp; 

temp = zeros(size(no_gsr,1),4,size(no_gsr,3),size(no_gsr,4),size(no_gsr,5));
temp(:,1,:,:,:) = squeeze(mean(no_gsr(:,1:2,:,:,:),2)); 
temp(:,2,:,:,:) = squeeze(mean(no_gsr(:,3:4,:,:,:),2)); 
temp(:,3,:,:,:) = squeeze(mean(no_gsr(:,5,:,:,:),2)); 
temp(:,4,:,:,:) = squeeze(mean(no_gsr(:,6,:,:,:),2)); 
no_gsr = temp; 

hrf = spm_hrf(0.693); 
zhrf = zeros(1,51); zhrf(26:end) = hrf(1:26) ;

state_colors = {[126,237,148]/255,[26,44,121]/255,[230,162,110]/255,[161,91,162]/255};
clear yesgsr_diffs nogsr_diffs
for i=1:5
   subplot(2,6,i) ;  
   plot(times,squeeze(mean(yes_gsr(:,3,3,i,:),1)),'LineWidth',2); hold on ; 
   plot(times,squeeze(mean(yes_gsr(:,4,3,i,:),1)),'LineWidth',2);  plot(times,zhrf,'k'); ylim([-.5,.25]);
   yesgsr_diffs(i,:) = mean(abs(squeeze(yes_gsr(:,4,3,i,:,1) - yes_gsr(:,3,3,i,:))),2); 
   mdiffs_yes(i) = mean(abs(squeeze(mean(yes_gsr(:,4,3,i,:),1)) - squeeze(mean(yes_gsr(:,3,3,i,:),1))));
   title(rsn_titles{i}); ylabel('correlation'); xlabel('time(s)'); xlim([-12,12]);
end
for i=1:5
   subplot(2,6,i+6) ;  
   plot(times,squeeze(mean(no_gsr(:,3,3,i,:),1)),'LineWidth',2); hold on ; 
   plot(times,squeeze(mean(no_gsr(:,4,3,i,:),1)),'LineWidth',2);  plot(times,zhrf,'k');   ylim([-.5,.25]);
   nogsr_diffs(i,:) = mean(abs(squeeze(no_gsr(:,4,3,i,:,1) - no_gsr(:,3,3,i,:))),2); 
   mdiffs_no(i) = mean(abs(squeeze(mean(no_gsr(:,4,3,i,:),1)) - squeeze(mean(no_gsr(:,3,3,i,:),1))));
   title(rsn_titles{i}); ylabel('correlation'); xlabel('time(s)'); xlim([-12,12]);
end
subplot(2,6,12); plot(1) ; hold on ; plot(2) ; plot(3,'k'); legend({'movie (eyes-open)','rest (eyes-closed)','canonical HRF'});

figure
subplot(1,2,1); 
b=bar(1,mean(mdiffs_no)); hold on ;
b.FaceColor = [1,0,0];
errorbar(1,mean(mdiffs_no),std(mdiffs_no)/sqrt(5),'k.');
b=bar(2,mean(mdiffs_yes)); hold on ;
errorbar(2,mean(mdiffs_yes),std(mdiffs_yes)/sqrt(5),'k.');
b.FaceColor = [0,1,0];
set(gca,'XTick',1:2,'XTickLabel',{'pre GSR','post GSR'});
[h,p,ci,stats] = ttest(mdiffs_no,mdiffs_yes);
title(['HRF difference ',format_p(p)]); ylabel('correlation'); 
subplot(1,2,2);
plot(1,'r','LineWidth',4) ; hold on ; plot(2,'g','LineWidth',4); legend({'pre-GSR','post-GSR'});







