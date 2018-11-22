clear all; close all; 
cd E:\russ_orient_data\russ_orientation_2

scans = {'russ_orientation_1.xdf','russ_orientation_2.xdf','russ_orientation_3.xdf','russ_orientation_4.xdf','russ_orientation_5.xdf'}; 

for scan=1:length(scans)
    disp(scan); 
eye0 = 'Pupil Primitive Data - Eye 0';
eye1 = 'Pupil Primitive Data - Eye 1'; 

rot2 = load_xdf(scans{scan}); 
clear names; 
for i=1:length(rot2)
   names{i} = rot2{i}.info.name ; 
   if strcmpi('StimulusMarkersLab',rot2{i}.info.name)
      marker_index = i;  
   end
end
ind0 = find(strcmpi(names,eye0)); 
ind1 = find(strcmpi(names,eye1)); 

ts1 = rot2{marker_index};

pup0 = rot2{ind0}.time_series(1,:); 
stamps0 = rot2{ind0}.time_series(3,:);
pup1 = rot2{ind1}.time_series(1,:); 
stamps1 = rot2{ind1}.time_series(3,:); 

trigs = ts1.time_series; 
trigstamps = ts1.time_stamps; 
stimtrigs = {'S2','S4','S12','S14','S22','S24','S32','S34','S42','S44'};
clear stimlats; 
for i=1:length(stimtrigs) ; lat_inds = find(strcmpi(stimtrigs{i},trigs)); stimlats(i,:) = trigstamps(lat_inds);  end
pup0 = interpolate_blinks(pup0,15,15); 
pup1 = interpolate_blinks(pup1,15,15); 

clear epochs0 epochs1; 
for i=1:size(stimlats,1)
    for j=1:size(stimlats,2)
        % eye 0
        diff_ij = abs(stimlats(i,j) - stamps0); 
        ind_ij = find(diff_ij==min(diff_ij),1); 
        epochs0(i,j,:) = pup0(ind_ij-300:ind_ij+100*24);    
        
        % eye 1
        diff_ij = abs(stimlats(i,j) - stamps1); 
        ind_ij = find(diff_ij==min(diff_ij),1); 
        epochs1(i,j,:) = pup1(ind_ij-300:ind_ij+100*24);    
    end
end
epochs0 = (epochs0 - repmat(mean(epochs0(:,:,1:300),3),[1,1,size(epochs1,3)])) ./ repmat(mean(epochs0(:,:,1:300),3),[1,1,size(epochs1,3)]);
epochs1 = (epochs1 - repmat(mean(epochs1(:,:,1:300),3),[1,1,size(epochs1,3)])) ./ repmat(mean(epochs1(:,:,1:300),3),[1,1,size(epochs1,3)]);

allepochs0(scan,:,:,:) = epochs0;
allepochs1(scan,:,:,:) = epochs1; 

end
contrast_colors = {[0.1,0.1,0],[0.2,0.2,0],[0.3,0.3,0],[0.4,0.4,0],[0.5,0.5,0]}; 
contrast_colors2 = {[0.1,0.1,0],[0.1,0.1,0],[0.2,0.2,0],[0.2,0.2,0],[0.3,0.3,0],[0.3,0.3,0],[0.4,0.4,0],[0.4,0.4,0],[0.5,0.5,0],[0.5,0.5,0]}; 
oblique_color = 'r'; vertical_color = 'c'; horizontal_color = 'g'; cardinal_color = 'm'; 

bothepochs(1,:,:,:,:) = allepochs0;
bothepochs(2,:,:,:,:) = allepochs1; 

rescounter = 1; 
resepochs = zeros(size(bothepochs,1),size(bothepochs,2)*size(bothepochs,4),size(bothepochs,3),size(bothepochs,5)); 
for i=1:size(bothepochs,3)
    rescounter = 1; 
    for j=1:size(bothepochs,2)
        for k=1:size(bothepochs,4)
            resepochs(:,rescounter,i,:) = squeeze(bothepochs(:,j,i,k,:)); 
            rescounter = rescounter + 1; 
        end
    end
end
times = (-300:2400)/100 ; 
for i=1:10
   %shadedErrorBar([],squeeze(mean(resepochs(1,:,i,:),2)),squeeze(std(resepochs(1,:,i,:),0,2))/sqrt(15));  hold on ; 
   plot(times,squeeze(mean(resepochs(1,:,i,:),2)),'Color',contrast_colors2{i},'LineWidth',3); hold on; 
end
vline([0,18],'r'); xlim([-3,24]); 
legend({'1','2','3','4','5','6','7','8','9','10'}); hline(0,'k'); 


s2_angles = [90:-1:0,359:-1:45];
s4_angles = [45:1:359,0:1:90];
time_inds = find(times > 0 & times < 18);    
gtime_inds = find(times>6 & times<18);    

tbersp = resepochs(:,:,:,time_inds); 
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
res_angles = mod(res_angles,180);

tbersp = resepochs(:,:,:,gtime_inds); 
res_angles = res_angles(:,(length(time_inds)-length(gtime_inds)+1):end); 

tbersp_amp = zeros(size(tbersp,1),size(tbersp,2),size(tbersp,3),180);    
for i=0:179
    for j=1:10
        bersp_ind = find(i==res_angles(j,:)); 
        tbersp_amp(:,:,j,i+1) = squeeze(mean(tbersp(:,:,j,bersp_ind),4)); 
    end
end

mbersp(1,:,:) = squeeze(mean(tbersp_amp(1,:,1:2,:),3)); 
mbersp(2,:,:) = squeeze(mean(tbersp_amp(1,:,3:4,:),3)); 
mbersp(3,:,:) = squeeze(mean(tbersp_amp(1,:,5:6,:),3)); 
mbersp(4,:,:) = squeeze(mean(tbersp_amp(1,:,7:8,:),3)); 
mbersp(5,:,:) = squeeze(mean(tbersp_amp(1,:,9:10,:),3)); 

for i=1:length(contrast_colors)
   shadedErrorBar(1:180,squeeze(mean(mbersp(i,:,:),2)),squeeze(std(mbersp(i,:,:),0,2))/sqrt(15),{'Color',contrast_colors{i}});  hold on ;
    
end


% plot bar charts for contrast 
oblique = [40:50,130:140];
vertical = 85:95; 
horizontal = [1:5,175:180];
cardinal = [85:95,1:5,175:180]; 
bold_contrast_labels = {'100','5'};

eyehigh = squeeze(mbersp(1,:,:))'; 
eyelow = squeeze(mbersp(5,:,:))'; 

subplot(4,3,1); 
for i=1:2:10
   shadedErrorBar(times,squeeze(mean(resepochs(1,:,i,:),2)),squeeze(std(resepochs(1,:,i,:),0,2))/sqrt(15),{'Color',contrast_colors2{i}});  hold on ; 
  % plot(times,squeeze(mean(resepochs(1,:,i,:),2)),'Color',contrast_colors2{i},'LineWidth',3); hold on; 
end
vline([0,18],'r'); xlim([-3,24]); ylim([-.25,.1]); hline(0,'k'); ylabel('%change'); xlabel('time(s)'); title('pupil diameter'); 
subplot(4,3,3) ; for i=1:5 ; plot(1,'Color',contrast_colors{i},'LineWidth',5) ; hold on ; end ; legend({'100% contrast','50% contrast','25% contrast','15% contrast','5% contrast'});

subplot(4,6,7);
b1 = bar(1,mean(mean(eyehigh(oblique,:)))); hold on; errorbar(1,mean(mean(eyehigh(oblique,:))),std(mean(eyehigh(oblique,:),1),0,2)/sqrt(14),'k.'); 
b1.FaceColor = oblique_color; 
b2 = bar(2,mean(mean(eyehigh(cardinal,:)))); errorbar(2,mean(mean(eyehigh(cardinal,:))),std(mean(eyehigh(cardinal,:),1),0,2)/sqrt(14),'k.'); 
b2.FaceColor = cardinal_color;
b3 = bar(4,mean(mean(eyelow(oblique,:)))); errorbar(4,mean(mean(eyelow(oblique,:))),std(mean(eyelow(oblique,:),1),0,2)/sqrt(14),'k.'); 
b3.FaceColor = oblique_color;
b4 = bar(5,mean(mean(eyelow(cardinal,:)))); errorbar(5,mean(mean(eyelow(cardinal,:))),std(mean(eyelow(cardinal,:),1),0,2)/sqrt(14),'k.'); 
b4.FaceColor = cardinal_color;
ymin = double(min([b1.YData, b2.YData, b3.YData, b4.YData])); ymax = double(max([b1.YData, b2.YData, b3.YData, b4.YData])); 
[~,p,~,~] = ttest(mean(eyehigh(oblique,:)),mean(eyehigh(cardinal,:))); if p<0.05 ; formp = ['*',format_p(p)]; else formp = format_p(p); end ;text(1.5,0.05,formp);
[~,p,~,~] = ttest(mean(eyelow(oblique,:)),mean(eyelow(cardinal,:))); if p<0.05 ; formp = ['*',format_p(p)]; else formp = format_p(p); end ; text(4.5,0.05,formp);
set(gca,'XTick',[1,4],'XTickLabel',bold_contrast_labels); xlabel('contrast(%)'); ylabel('%change'); title('oblique vs cardinal'); ylim([-0.15,0.075]);

subplot(4,6,8); 
b1 = bar(1,mean(mean(eyehigh(oblique,:)))); hold on; errorbar(1,mean(mean(eyehigh(oblique,:))),std(mean(eyehigh(oblique,:),1),0,2)/sqrt(14),'k.'); 
b1.FaceColor = oblique_color; 
b2 = bar(2,mean(mean(eyehigh(vertical,:)))); errorbar(2,mean(mean(eyehigh(vertical,:))),std(mean(eyehigh(vertical,:),1),0,2)/sqrt(14),'k.'); 
b2.FaceColor = vertical_color; 
b3 = bar(4,mean(mean(eyelow(oblique,:)))); errorbar(4,mean(mean(eyelow(oblique,:))),std(mean(eyelow(oblique,:),1),0,2)/sqrt(14),'k.'); 
b3.FaceColor = oblique_color;
b4 = bar(5,mean(mean(eyelow(vertical,:)))); errorbar(5,mean(mean(eyelow(vertical,:))),std(mean(eyelow(vertical,:),1),0,2)/sqrt(14),'k.'); 
b4.FaceColor = vertical_color; 
ymin = double(min([b1.YData, b2.YData, b3.YData, b4.YData])); ymax = double(max([b1.YData, b2.YData, b3.YData, b4.YData])); 
[~,p,~,~] = ttest(mean(eyehigh(oblique,:)),mean(eyehigh(vertical,:))); if p<0.05 ; formp = ['*',format_p(p)]; else formp = format_p(p); end ;text(1.5,0.05,formp);
[~,p,~,~] = ttest(mean(eyelow(oblique,:)),mean(eyelow(vertical,:))); if p<0.05 ; formp = ['*',format_p(p)]; else formp = format_p(p); end ; text(4.5,0.05,formp);
set(gca,'XTick',[1,4],'XTickLabel',bold_contrast_labels); xlabel('contrast(%)'); ylabel('%change'); title('oblique vs vertical');ylim([-0.15,0.075]);

subplot(4,6,9); 
b1 = bar(1,mean(mean(eyehigh(vertical,:)))); hold on; errorbar(1,mean(mean(eyehigh(vertical,:))),std(mean(eyehigh(vertical,:),1),0,2)/sqrt(14),'k.'); 
b1.FaceColor = vertical_color; 
b2 = bar(2,mean(mean(eyehigh(horizontal,:)))); errorbar(2,mean(mean(eyehigh(horizontal,:))),std(mean(eyehigh(horizontal,:),1),0,2)/sqrt(14),'k.'); 
b2.FaceColor = horizontal_color; 
b3 = bar(4,mean(mean(eyelow(vertical,:)))); errorbar(4,mean(mean(eyelow(vertical,:))),std(mean(eyelow(vertical,:),1),0,2)/sqrt(14),'k.'); 
b3.FaceColor = vertical_color;
b4 = bar(5,mean(mean(eyelow(horizontal,:)))); errorbar(5,mean(mean(eyelow(horizontal,:))),std(mean(eyelow(horizontal,:),1),0,2)/sqrt(14),'k.'); 
b4.FaceColor = horizontal_color; 
ymin = double(min([b1.YData, b2.YData, b3.YData, b4.YData])); ymax = double(max([b1.YData, b2.YData, b3.YData, b4.YData])); 
[~,p,~,~] = ttest(mean(eyehigh(vertical,:)),mean(eyehigh(horizontal,:))); if p<0.05 ; formp = ['*',format_p(p)]; else formp = format_p(p); end ;text(1.5,0.05,formp);
[~,p,~,~] = ttest(mean(eyelow(vertical,:)),mean(eyelow(horizontal,:))); if p<0.05 ; formp = ['*',format_p(p)]; else formp = format_p(p); end ; text(4.5,0.05,formp);
set(gca,'XTick',[1,4],'XTickLabel',bold_contrast_labels); xlabel('contrast(%)'); ylabel('%change'); title('vertical vs horizontal');ylim([-0.15,0.075]);

subplot(4,6,13) ; plot(1,oblique_color,'LineWidth',5); hold on ; plot(1,cardinal_color,'LineWidth',5); legend({'oblique','cardinal'});
subplot(4,6,14) ; plot(1,oblique_color,'LineWidth',5); hold on ; plot(1,vertical_color,'LineWidth',5); legend({'oblique','vertical'});
subplot(4,6,15) ; plot(1,vertical_color,'LineWidth',5); hold on ; plot(1,horizontal_color,'LineWidth',5); legend({'vertical','horizontal'});

set(gcf,'Position',[100 100 1700 900])











