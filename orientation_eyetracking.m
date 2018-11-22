clear all ; close all; 
subs = {'russell','greg'};
pooldata = zeros(10,39,5632); 
for sb=1:length(subs)
cd(['E:\orientation_eyetracking\',subs{sb}]);

trigs = {'S2.','S4.','S12.','S14.','S22.','S24.','S32.','S34.','S42.','S44.'};
clear alldata; 
for tr=1:length(trigs)
    
    name = dir(['*eye0*',trigs{tr},'*']);
    data = load(name.name); 
    
    data = eegfiltfft(data,256,0,4); 
    
    alldata(tr,:,:) = data; 
    

end
CC = {[194,82,60]/255,[230,142,28]/255,[123,237,0]/255,[30,144,148]/255,[11,44,122]/255}; 
CONTRAST_COLORS = {CC{1},CC{1},CC{2},CC{2},CC{3},CC{3},CC{4},CC{4},CC{5},CC{5}}; 
cols1 = [1:90,90:-1:1];
cols2 = [1:45,45:-1:1,1:45,45:-1:1]; 
cols3 = [90:-1:1,1:90];
colsrgb(1,:,1) = mat2gray(cols1); colsrgb(1,:,2) = mat2gray(cols3); colsrgb(1,:,3) = mat2gray(cols2); 
oblique_color = squeeze(colsrgb(1,45,:));
vertical_color = squeeze(colsrgb(1,90,:)); 
horizontal_color = squeeze(colsrgb(1,1,:)); 
cardinal_color = squeeze(colsrgb(1,1,:))/2 + squeeze(colsrgb(1,90,:))/2; 

alldata = (alldata - repmat(mean(alldata(:,:,1:512),3),[1,1,size(alldata,3)])) ./ repmat(mean(alldata(:,:,1:512),3),[1,1,size(alldata,3)]); 
size(alldata)
if sb==1 ; pooldata(:,1:30,:) = alldata; else pooldata(:,31:39,:) = alldata; end

end
alldata = pooldata; 

labels = {'100% (CW)','100% (CCW)','50% (CW)','50% (CCW)','25% (CW)','25% (CCW)','15% (CW)','15% (CCW)','5% (CW)','5% (CCW)'};

subplot(4,3,1); 
times = [(-256*2):(256*20)-1]/256;
for i=1:10
   plot(times,squeeze(mean(alldata(i,:,:),2)),'Color',CONTRAST_COLORS{i},'LineWidth',3); hold on; 
end
legend(labels); hline(0,'k'); vline([0,18]); xlabel('time(s)'); ylabel('pupil diameter %change'); xlim([-2,20]); ylim([-.15,.05]);

s2_angles = [90:-1:0,359:-1:45];
s4_angles = [45:1:359,0:1:90];
time_inds = find(times > 0 & times < 18);    
gtime_inds = find(times>6 & times<18);    

tbersp = alldata(:,:,time_inds); 
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

tbersp = alldata(:,:,gtime_inds); 
res_angles = res_angles(:,(length(time_inds)-length(gtime_inds)+1):end); 
clear tbersp_amp
for i=0:179
    for j=1:10
        bersp_ind = find(i==res_angles(j,:)); 
        tbersp_amp(:,j,i+1) = squeeze(mean(tbersp(j,:,bersp_ind),3)); 
    end
end

mbersp(1,:,:) = squeeze(mean(tbersp_amp(:,1:2,:),2)); 
mbersp(2,:,:) = squeeze(mean(tbersp_amp(:,3:4,:),2)); 
mbersp(3,:,:) = squeeze(mean(tbersp_amp(:,5:6,:),2)); 
mbersp(4,:,:) = squeeze(mean(tbersp_amp(:,7:8,:),2)); 
mbersp(5,:,:) = squeeze(mean(tbersp_amp(:,9:10,:),2)); 
%{
for i=1:length(CONTRAST_COLORS)
   shadedErrorBar(1:180,squeeze(mean(mbersp(i,:,:),2)),squeeze(std(mbersp(i,:,:),0,2))/sqrt(15),{'Color',CONTRAST_COLORS{i}});  hold on ;
    
end
%}

% plot bar charts for contrast 
oblique = [40:50,130:140];
vertical = 85:95; 
horizontal = [1:5,175:180];
cardinal = [85:95,1:5,175:180]; 
bold_contrast_labels = {'100','5'};

eyehigh = squeeze(mbersp(1,:,:))'; 
eyelow = squeeze(mbersp(5,:,:))'; 



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
set(gca,'XTick',[1,4],'XTickLabel',bold_contrast_labels); xlabel('contrast(%)'); ylabel('%change'); title('oblique vs cardinal'); ylim([-0.1,0.075]);

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
set(gca,'XTick',[1,4],'XTickLabel',bold_contrast_labels); xlabel('contrast(%)'); ylabel('%change'); title('oblique vs vertical');ylim([-0.1,0.075]);

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
set(gca,'XTick',[1,4],'XTickLabel',bold_contrast_labels); xlabel('contrast(%)'); ylabel('%change'); title('vertical vs horizontal');ylim([-0.1,0.075]);

subplot(4,6,13) ; plot(1,'Color',oblique_color,'LineWidth',5); hold on ; plot(1,'Color',cardinal_color,'LineWidth',5); legend({'oblique','cardinal'});
subplot(4,6,14) ; plot(1,'Color',oblique_color,'LineWidth',5); hold on ; plot(1,'Color',vertical_color,'LineWidth',5); legend({'oblique','vertical'});
subplot(4,6,15) ; plot(1,'Color',vertical_color,'LineWidth',5); hold on ; plot(1,'Color',horizontal_color,'LineWidth',5); legend({'vertical','horizontal'});

set(gcf,'Position',[100 100 1700 900])

