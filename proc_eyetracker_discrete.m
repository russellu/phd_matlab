clear all; close all; 
cd E:\russ_eye_date_07122018

scans = {'discrete_1.xdf','discrete_2.xdf','discrete_3.xdf'}; 

for scan=1:length(scans); 
eye0 = 'Pupil Primitive Data - Eye 0';
eye1 = 'Pupil Primitive Data - Eye 1'; 

rot2 = load_xdf(scans{scan}); 
clear names; 
for i=1:length(rot2)
   names{i} = rot2{i}.info.name ; 
end
ind0 = find(strcmpi(names,eye0)); 
ind1 = find(strcmpi(names,eye1)); 

ts1 = rot2{1};

pup0 = rot2{ind0}.time_series(1,:); 
stamps0 = rot2{ind0}.time_series(3,:);
pup1 = rot2{ind1}.time_series(1,:); 
stamps1 = rot2{ind1}.time_series(3,:); 

trigs = ts1.time_series; 
trigstamps = ts1.time_stamps; 
stimtrigs = {'S1','S2','S3','S4','S5','S14','S15','S16','S17','S19','S20','S31','S32','S33','S34','S35'};
for i=1:length(stimtrigs) ; lat_inds = find(strcmpi(stimtrigs{i},trigs)); stimlats(i,:) = trigstamps(lat_inds);  end

clear epochs0 epochs1; 
for i=1:size(stimlats,1)
    for j=1:size(stimlats,2)
        % eye 0
        diff_ij = abs(stimlats(i,j) - stamps0); 
        ind_ij = find(diff_ij==min(diff_ij),1); 
        epochs0(i,j,:) = pup0(ind_ij-90:ind_ij+60*4);    
        
        % eye 1
        diff_ij = abs(stimlats(i,j) - stamps1); 
        ind_ij = find(diff_ij==min(diff_ij),1); 
        epochs1(i,j,:) = pup1(ind_ij-90:ind_ij+60*4);    
    end
end
epochs0 = (epochs0 - repmat(mean(epochs0(:,:,1:90),3),[1,1,331])) ./ repmat(mean(epochs0(:,:,1:90),3),[1,1,331]);
epochs1 = (epochs1 - repmat(mean(epochs1(:,:,1:90),3),[1,1,331])) ./ repmat(mean(epochs1(:,:,1:90),3),[1,1,331]);

allepochs0(scan,:,:,:) = epochs0; allepochs1(scan,:,:,:) = epochs1; 

end

times=(-90:60*4)/60; 

for i=1:3 ; for j=1:16; for k=1:10 ; plot(squeeze(allepochs0(i,j,k,:))) ; hold on; end ; end ; end
% combine eyes into a single vector:
bothepochs(1,:,:,:,:) = allepochs0;
bothepochs(2,:,:,:,:) = allepochs1; 
% combine all trials:
resboth = zeros(size(bothepochs,1),size(bothepochs,3),size(bothepochs,2)*size(bothepochs,3),size(bothepochs,5)); 
for i=1:size(bothepochs,1)
    for j=1:size(bothepochs,3)
        resboth(i,j,:,:) = reshape(squeeze(bothepochs(i,:,j,:,:)),[size(bothepochs,2)*size(bothepochs,4),331]);
    end
end

resepochs = reshape(resboth,[numel(resboth(:,:,:,1)),size(resboth,4)]);
zepochs = (zscore(sum(abs(diff(resepochs,1,2)),2))) < 0;
zepochs = reshape(zepochs,[size(resboth(:,:,:,1))]); 
    
titles = {'full','left','right','top','bottom','bottom right','bottom left','top left','top right','periphery','fovea','1%','5%','15%','30%','100%'};

subplot(2,2,1);
colors = {[1,0,1],[.7,.7,0],[.4,.8,.8],[0,1,0],[0,0,0]};
for i=12:16  
    plot(times,squeeze(mean(resboth(2,i,zepochs(2,i,:),:),3)),'Color',colors{i-11}); hold on ; 
end
legend(titles(12:16)); 
for i=12:16
    shadedErrorBar(times,squeeze(mean(resboth(2,i,zepochs(2,i,:),:),3)),squeeze(std(resboth(2,i,zepochs(2,i,:),:),0,3))/sqrt(30),{'Color',colors{i-11}}); hold on ; 
    hline(0,'k') ; vline(0,'k') ; vline(1.75,'k'); 
    ylim([-.3,.05]);
end
xlabel('time(s)'); ylabel('diameter %change'); title('contrast'); 

subplot(2,2,2);
colors = {[1,0,1],[.7,.7,0],[.4,.8,.8],[0,1,0],[0,0,0]};
for i=1:5  
    plot(times,squeeze(mean(resboth(2,i,zepochs(2,i,:),:),3)),'Color',colors{i}); hold on ; 
end
legend(titles(1:5)); 
for i=1:5
    shadedErrorBar(times,squeeze(mean(resboth(2,i,zepochs(2,i,:),:),3)),squeeze(std(resboth(2,i,zepochs(2,i,:),:),0,3))/sqrt(30),{'Color',colors{i}}); hold on ; 
    hline(0,'k') ; vline(0,'k') ; vline(1.75,'k'); 
    ylim([-.3,.05]);
end
xlabel('time(s)'); ylabel('diameter %change'); title('hemis + full'); 


subplot(2,2,3);
colors = {[1,0,1],[.7,.7,0],[.4,.8,.8],[0,1,0],[0,0,0]};
for i=6:9
    plot(times,squeeze(mean(resboth(2,i,zepochs(2,i,:),:),3)),'Color',colors{i-5}); hold on ; 
end
legend(titles(6:9)); 
for i=6:9
    shadedErrorBar(times,squeeze(mean(resboth(2,i,zepochs(2,i,:),:),3)),squeeze(std(resboth(2,i,zepochs(2,i,:),:),0,3))/sqrt(30),{'Color',colors{i-5}}); hold on ; 
    hline(0,'k') ; vline(0,'k') ; vline(1.75,'k'); 
    ylim([-.1,.05]);
end
xlabel('time(s)'); ylabel('diameter %change'); title('quadrants'); 

subplot(2,2,4);
colors = {[1,0,1],[.7,.7,0],[.4,.8,.8],[0,1,0],[0,0,0]};
for i=10:11
    plot(times,squeeze(mean(resboth(2,i,zepochs(2,i,:),:),3)),'Color',colors{i-9}); hold on ; 
end
legend(titles(10:11)); 
for i=10:11
    shadedErrorBar(times,squeeze(mean(resboth(2,i,zepochs(2,i,:),:),3)),squeeze(std(resboth(2,i,zepochs(2,i,:),:),0,3))/sqrt(30),{'Color',colors{i-9}}); hold on ; 
    hline(0,'k') ; vline(0,'k') ; vline(1.75,'k'); 
    ylim([-.1,.05]);
end
xlabel('time(s)'); ylabel('diameter %change'); title('periphery + fovea'); 



%{
subplot(2,2,1);
plot(squeeze(mean(epochs(1:5,:,:),2))') ; vline(90,'k'); hline(0,'k'); legend(titles(1:5)); title('full + hemis'); 
subplot(2,2,2);
plot(squeeze(mean(epochs(6:9,:,:),2))') ; vline(90,'k'); hline(0,'k'); legend(titles(6:9)); title('quadrants');
subplot(2,2,3);
plot(squeeze(mean(epochs(10:11,:,:),2))') ; vline(90,'k'); hline(0,'k'); legend(titles(10:11)); title('periphery vs fovea');
subplot(2,2,4);
plot(squeeze(mean(epochs(12:16,:,:),2))') ; vline(90,'k'); hline(0,'k'); legend(titles(12:16)); title('contrast');
%}