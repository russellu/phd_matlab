clear all ; close all ;
rsns = {'anterior_visual.nii','dmn.nii','higher_visual.nii','primary_visual.nii','motor.nii'};
rsn_titles = {'Anterior Visual','Default Mode','Higher-order Visual','Primary Visual','Motor'};
cd C:\shared\epireg\rsns\networks;

for r=1:length(rsns)
   rsn = load_untouch_nii(rsns{r}); 
   [sv,si] = sort(rsn.img(:),'descend'); 
   rsnvox(r,:) = si(1:1000);    
end

zimg = zeros(size(rsn.img)); 

cd E:\sim\saved

nogsr = load('mc_no_gsr'); nogsr = nogsr.mc_no_gsr; 
temp = zeros(size(nogsr,1),size(nogsr,2),size(nogsr,3),6,4,7); 
temp(:,:,:,:,1,:) = mean(nogsr(:,:,:,:,1:2,:),5); 
temp(:,:,:,:,2,:) = mean(nogsr(:,:,:,:,3:4,:),5); 
temp(:,:,:,:,3,:) = nogsr(:,:,:,:,5,:); 
temp(:,:,:,:,4,:) = nogsr(:,:,:,:,6,:); 
nogsr = temp; 

yesgsr = load('mc_yes_gsr'); yesgsr = yesgsr.mc_yes_gsr; 
temp = zeros(size(yesgsr,1),size(yesgsr,2),size(yesgsr,3),6,4,7); 
temp(:,:,:,:,1,:) = mean(yesgsr(:,:,:,:,1:2,:),5); 
temp(:,:,:,:,2,:) = mean(yesgsr(:,:,:,:,3:4,:),5); 
temp(:,:,:,:,3,:) = yesgsr(:,:,:,:,5,:); 
temp(:,:,:,:,4,:) = yesgsr(:,:,:,:,6,:); 
yesgsr = temp; 

slices = [5,8,10,12,14,16,20]; 
cd E:\sim\fmri\alex ; meanatlas = load_untouch_nii('mean.nii.gz'); 

thresh = 0.1; 
for sl=1:length(slices)
    subplottight(1,7,sl); 
    corrimg = squeeze(mean(mean(yesgsr(:,:,slices(sl),5,2,:),3),6)); corrimg(1,1) = -thresh; corrimg(1,2) = thresh; corrimg(corrimg>thresh) = thresh; corrimg(corrimg<-thresh) = -thresh; 
    corrimg(isnan(corrimg)) = 0; 
    plotoverlayIntensity2D(meanatlas.img(:,:,slices(sl)),mat2gray(abs(corrimg)),corrimg,270);    
end

res_nogsr = reshape(nogsr,[64*64*33,6,4,7]);
res_yesgsr = reshape(yesgsr,[64*64*33,6,4,7]); 
clear rs_nogsr rs_yesgsr; 
for rs=1:length(rsns)
   rs_nogsr(rs,:,:,:) = squeeze(mean(res_nogsr(rsnvox(rs),:,:,:),1));  
   rs_yesgsr(rs,:,:,:) = squeeze(mean(res_yesgsr(rsnvox(rs),:,:,:),1)); 
end

for i=1:5
subplot(1,6,i); 
rsn_ind = i; freq_ind = 5; state_ind = 2;  
b = bar(1,squeeze(mean(rs_nogsr(rsn_ind,freq_ind,state_ind,:),4))) ; hold on ; 
b.FaceColor = [1,0,0]; 
errorbar(1,squeeze(mean(rs_nogsr(rsn_ind,freq_ind,state_ind,:),4)),squeeze(std(rs_nogsr(rsn_ind,freq_ind,state_ind,:),0,4))/sqrt(7),'k.'); 
b = bar(2,squeeze(mean(rs_yesgsr(rsn_ind,freq_ind,state_ind,:),4))) ; 
b.FaceColor = [0,1,0]; 
errorbar(2,squeeze(mean(rs_yesgsr(rsn_ind,freq_ind,state_ind,:),4)),squeeze(std(rs_yesgsr(rsn_ind,freq_ind,state_ind,:),0,4))/sqrt(7),'k.'); 
[h,p,ci,stats] = ttest(squeeze(rs_nogsr(rsn_ind,freq_ind,state_ind,:)),squeeze(rs_yesgsr(rsn_ind,freq_ind,state_ind,:))); 
title([rsn_titles{rsn_ind},' ',format_p(p)]); ylabel('correlation (rho)'); set(gca,'XTick',[1,2],'XTickLabel',{'pre GSR','post GSR'});
end

subplot(1,6,6); plot(1,'r','LineWidth',4); hold on ; plot(2,'g','LineWidth',4); legend({'pre GSR','post GSR'})


state_colors = {[126,237,148]/255,[26,44,121]/255,[230,162,110]/255,[161,91,162]/255};
for rsn=1:5
    subplot(2,6,rsn);
    for i=1:4
       freq_ind = 5; rsn_ind = rsn; 
       b = bar(i,squeeze(mean(rs_nogsr(rsn_ind,freq_ind,i,:),4))); hold on; 
       b.FaceColor = state_colors{i}; 
       errorbar(i,squeeze(mean(rs_nogsr(rsn_ind,freq_ind,i,:),4)),squeeze(std(rs_nogsr(rsn_ind,freq_ind,i,:),0,4))/sqrt(7),'k.')
    end
    p = anova1(squeeze(rs_nogsr(rsn_ind,freq_ind,:,:))',[],'off');
    title([rsn_titles{rsn},' ',format_p(p)]);
    set(gca,'XTick',1:4,'XTickLabel',[]) ; ylabel('correlation (rho)'); xlabel('RSN');  
end
subplot(2,6,12); plot(1,'Color',state_colors{1},'LineWidth',3); hold on; plot(1,'Color',state_colors{2},'LineWidth',3);plot(1,'Color',state_colors{3},'LineWidth',3);plot(1,'Color',state_colors{4},'LineWidth',3);
legend({'retinotopy','event-related','movie','rest'});
for rsn=1:5
    subplot(2,6,rsn+6);
    for i=1:4
       freq_ind = 5; rsn_ind = rsn; 
       b = bar(i,squeeze(mean(rs_yesgsr(rsn_ind,freq_ind,i,:),4))); hold on; 
       b.FaceColor = state_colors{i}; 
       errorbar(i,squeeze(mean(rs_yesgsr(rsn_ind,freq_ind,i,:),4)),squeeze(std(rs_yesgsr(rsn_ind,freq_ind,i,:),0,4))/sqrt(7),'k.')
    end
    p = anova1(squeeze(rs_yesgsr(rsn_ind,freq_ind,:,:))',[],'off');
    title([rsn_titles{rsn},' ',format_p(p)]);
    set(gca,'XTick',1:4,'XTickLabel',[]) ; ylabel('correlation (rho)'); xlabel('RSN');  
end
subplot(2,6,12); plot(1,'Color',state_colors{1},'LineWidth',3); hold on; plot(1,'Color',state_colors{2},'LineWidth',3);plot(1,'Color',state_colors{3},'LineWidth',3);plot(1,'Color',state_colors{4},'LineWidth',3);
legend({'retinotopy','event-related','movie','rest'});












