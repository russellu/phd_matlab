clear all ; close all; 
cd e:/sim/saved; 
% xcorrs(sub,scan,freq,rsn,:)

rsn_titles = {'Anterior Visual','Default Mode','Higher-order Visual','Primary Visual','Motor'};

nogsr = load('sortpower_gsr') ; nogsr = nogsr.sortpower_gsr; 

nogsr(isnan(nogsr)) = 0; 
temp = zeros(size(nogsr,1),4,size(nogsr,3),size(nogsr,4),size(nogsr,5)); 
temp(:,1,:,:,:) = mean(nogsr(:,1:2,:,:,:),2); 
temp(:,2,:,:,:) = mean(nogsr(:,3:4,:,:,:),2); 
temp(:,3,:,:,:) = mean(nogsr(:,5,:,:,:),2); 
temp(:,4,:,:,:) = mean(nogsr(:,6,:,:,:),2); 
nogsr = temp; 

freq_ind = 3; task_ind =4; 
topinds = 1:40; botinds = 400-40:400; 
for i=1:5
    subplot(1,6,i)
    h = bar(1,mean(mean(nogsr(:,task_ind,freq_ind,i,topinds),1),5)); hold on; 
    h.FaceColor = [1,0.5,0.5];
    errorbar(1,mean(mean(nogsr(:,task_ind,freq_ind,i,topinds),1),5),std(mean(nogsr(:,task_ind,freq_ind,i,topinds),5),0,1)/sqrt(7),'k.'); 
    
    h = bar(2,mean(mean(nogsr(:,task_ind,freq_ind,i,botinds),1),5));
    h.FaceColor = [0.5,1,0.5];
    errorbar(2,mean(mean(nogsr(:,task_ind,freq_ind,i,botinds),1),5),std(mean(nogsr(:,task_ind,freq_ind,i,botinds),5),0,1)/sqrt(7),'k.'); 
    
    [h,p,ci,stats] = ttest(squeeze(mean(nogsr(:,task_ind,freq_ind,i,topinds),5)),squeeze(mean(nogsr(:,task_ind,freq_ind,i,botinds),5)));
    title([rsn_titles{i},' ',format_p(p)]);
end


straightcorrs = load('straightcorrs_gsr') ; straightcorrs = straightcorrs.straightcorrs_gsr; 
temp = zeros(size(straightcorrs,1),4,size(straightcorrs,3),size(straightcorrs,4)); 
temp(:,1,:,:) = mean(straightcorrs(:,1:2,:,:),2); 
temp(:,2,:,:) = mean(straightcorrs(:,3:4,:,:),2); 
temp(:,3,:,:) = mean(straightcorrs(:,5,:,:),2); 
temp(:,4,:,:) = mean(straightcorrs(:,6,:,:),2); 
straightcorrs = temp ;


bar(squeeze(mean(straightcorrs(:,1,3,:),1))); 












