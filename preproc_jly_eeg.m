clear all ; close all
subs = {'reihaneh','russell','camille','ariane','annfred','claudie','gabriel','krisen','Jessily'};
bades = {[56],[16,19],[45],[6],[],[16],[32],[16,19,57],[42,45,12]};
bestcomps ={[],[20],[9],[6],[3],[4],[7],[15],[14]}; 
flips = [0,0,0,0,0,0,0,1,1];

cd e:/saved ; discrete_template = load('discrete_template'); discrete_template = discrete_template.discrete_template;
dt = squeeze(mean(mean(discrete_template,1),2)); 
gamma_dt = squeeze(mean(dt(15:40,:),1)); alpha_dt = squeeze(mean(dt(4:12,:),1)); 

for sb=1:length(subs) 
   
cd(['E:\jly_russ\',subs{sb}]);

%{
vhdrs = dir('*vhdr'); 

for filename = 1:length(vhdrs)
eeg = pop_loadbv('.',vhdrs(filename).name);  % load the raw file
eeg = pop_chanedit(eeg,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
eeg = pop_resample(eeg,250); 
eeg.data = eeg.data - eegfiltfft(eeg.data,eeg.srate,59,61); % new data = old data minus 60Hz
eeg.data = eeg.data - eegfiltfft(eeg.data,eeg.srate,84,86); 
if filename==1; merged = eeg; else; merged = pop_mergeset(merged,eeg); end 
end

if flips(sb) == 1 ; temp = merged.data(1:32,:); merged.data(1:32,:) = merged.data(33:64,:); merged.data(33:64,:) = temp ;end
pop_saveset(merged,'mergefilt.set'); 

eeg = pop_loadset('mergefilt.set'); 
figure,bar(zscore(sum(abs(diff(eeg.data,1,2)),2)))
merged = pop_interp(eeg,bades{sb},'spherical');
pop_saveset(eeg,'mergefilt.set'); 
%}
eeg = pop_loadset('mergefilt.set'); 

%{
eeg = pop_loadset('mergefilt.set'); 
highpass_filtered = eegfiltfft(eeg.data,eeg.srate,1,120); % remove low frequencies
[weights, sphere] = runica(highpass_filtered,'maxsteps',128); % run ICA
winv = pinv(weights*sphere);  % inverse weight matrix (to go back to electrodes from components)
acts = weights*sphere*eeg.data; % ICA component time courses (activations)
figure,for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),eeg.chanlocs);  end ;suptitle(subs{sb});
saveweights{1} = weights; saveweights{2} = sphere; 
save('saveweights','saveweights'); 
%}


saveweights = load('saveweights'); saveweights = saveweights.saveweights;
weights = saveweights{1}; sphere = saveweights{2};
winv = pinv(weights*sphere); 
acts = weights*sphere*eeg.data; 
%figure,for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),merged.chanlocs);  end ;suptitle(subs{sb});

neweeg = eeg; % copy eeg structure
neweeg.data = acts; % replace channel data with activations



%{
trigs = {'S  1','S  2','S  3','S  4','S  5','S 14','S 15','S 16','S 17','S 19','S 20','S 31','S 32','S 33','S 34','S 35'}; 
clear allersp 
epochs = pop_epoch(neweeg,trigs,[-.5,2.25]); % epoch data around stimulus markers
for i=1:64 % run time-frequency analysis on all components (64)
    disp(i); jcount = 1; 
    for j=1:2:size(epochs.data,3)
        [allersp(i,jcount,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochs.data(i,:,j)),epochs.pnts,[epochs.xmin,epochs.xmax],epochs.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,100],'nfreqs',50,'winsize',64,'baseline',NaN,'verbose','off','timesout',100) ; 
            jcount = jcount + 1; 
    end
end
clear alpha_corrs gamma_corrs;
bersp = allersp - repmat(mean(allersp(:,:,:,times<0),4),[1,1,1,100]); 
gamma_bersp = squeeze(mean(bersp(:,:,15:40,:),3));
alpha_bersp = squeeze(mean(bersp(:,:,4:12,:),3)); 
for i=1:64 ; alpha_corrs(i,:) = corr(squeeze(alpha_bersp(i,:,:))',alpha_dt') ; end
for i=1:64 ; gamma_corrs(i,:) = corr(squeeze(gamma_bersp(i,:,:))',gamma_dt') ; end

[sv,si] = sort(mean(alpha_corrs,2) + mean(gamma_corrs,2),'descend'); 
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(bersp(si(i),:,:,:),2)),[-2,2]) ; axis xy ; end
sorted_comp_indices = si ; save('sorted_comp_indices','sorted_comp_indices'); 

%}







si = load('sorted_comp_indices'); si =si.sorted_comp_indices; 
winv(:,si(6:end)) = 0; 
invacts = winv*acts; neweeg.data = invacts; 
trigs = {'S  1','S  2','S  3','S  4','S  5','S 14','S 15','S 16','S 17','S 19','S 20','S 31','S 32','S 33','S 34','S 35'}; 
gamma_invacts = eegfiltfft(invacts,eeg.srate,30,90); 
alpha_invacts = eegfiltfft(invacts,eeg.srate,8,25); 

neweeg.data = gamma_invacts; 
for tr=1:length(trigs)
    epochs = pop_epoch(neweeg,{trigs{tr}},[-.5,2.25]); % epoch data around stimulus markers
    epochs.data = abs(epochs.data); 
    gamma_elecs(sb,tr,:) = squeeze(mean(mean(epochs.data(:,epochs.times>0 & epochs.times<1500,:),2) - mean(epochs.data(:,epochs.times<0,:),2),3));    
end

neweeg.data = alpha_invacts; 
for tr=1:length(trigs)
    epochs = pop_epoch(neweeg,{trigs{tr}},[-.5,2.25]); % epoch data around stimulus markers
    epochs.data = abs(epochs.data); 
    alpha_elecs(sb,tr,:) = squeeze(mean(mean(epochs.data(:,epochs.times>0 & epochs.times<1500,:),2) - mean(epochs.data(:,epochs.times<0,:),2),3));    
end



trigs = {'S  1','S  2','S  3','S  4','S  5','S 14','S 15','S 16','S 17','S 19','S 20','S 31','S 32','S 33','S 34','S 35'}; 
clear allersp 
neweeg.data = acts; 
for tr=1:length(trigs)
epochs = pop_epoch(neweeg,{trigs{tr}},[-.85,2.85]); % epoch data around stimulus markers
badts = find(zscore(squeeze(mean(std(epochs.data(:,1:200,:),0,2),1))) > 1.5);
epochs.data(:,:,badts) = []; 
for i=1:5 % run time-frequency analysis on all components (64)
        [allersp(tr,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochs.data(si(i),:,:)),epochs.pnts,[epochs.xmin,epochs.xmax],epochs.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,100],'nfreqs',50,'winsize',64,'baseline',0,'verbose','off','timesout',300) ; 
end
end
figure,for i=1:5 ; subplot(1,5,i) ; imagesc(squeeze(mean(allersp(:,i,:,:),1)),[-2,2]) ; end ; suptitle(subs{sb}); 

mersp(sb,:,:,:) = squeeze(mean(allersp(:,1:3,:,:),2)) ;


%}
%comps = subcomps{sb}; 
%subersp(sb,:,:,:) = squeeze(mean(allersp(:,comps,:,:),2)); 
%figure,
%for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(allersp(:,i,:,:),1)),[-2,2]) ; axis xy  ;colormap jet; title(i); end ; suptitle(subs{sb}); 
%subersp(sb,:,:,:) = squeeze(allersp); 


%{
%}
end



for i=1:11 ; subplot(4,12,i+12) ; 
    imagesc(times,freqs,squeeze(mean(mersp(:,i,:,:),1)),[-3,3]) ; axis xy ;
    if i==1 ; xlabel('time(s)'); ylabel('frequency(hz)'); vline([0,1.75],'k'); colormap jet; else ; set(gca,'XTickLabel',[],'YTickLabel',[]) ; end; 
end
for i=1:11 ; subplot(4,12,i+24) ; topoplot(squeeze(mean(gamma_elecs(:,i,:),1)),eeg.chanlocs,'maplimits',[-.3,.3],'style','both','electrodes','off') ; end ; colormap jet;
for i=1:11 ; subplot(4,12,i+36) ; topoplot(squeeze(mean(alpha_elecs(:,i,:),1)),eeg.chanlocs,'maplimits',[-1.5,1.5],'style','both','electrodes','off') ; end ; colormap jet; 
subplot(4,12,24) ; imagesc([-3,3]) ; h=colorbar; title(h,'dB'); 
subplot(4,12,36) ; imagesc([-.3,.3]) ; h=colorbar; title(h,'\muV (gamma)'); 
subplot(4,12,48) ; imagesc([-1.5,1.5]) ; h=colorbar; title(h,'\muV (alpha/beta)'); 

figure, contrastinds = 16:-1:12; 
for i=1:5 ; subplot(4,12,i) ; 
    imagesc(times,freqs,squeeze(mean(mersp(:,contrastinds(i),:,:),1)),[-3,3]) ; axis xy ;
    if i==1 ; xlabel('time(s)'); ylabel('frequency(hz)'); vline([0,1.75],'k'); colormap jet; else ; set(gca,'XTickLabel',[],'YTickLabel',[]) ; end; 
end
for i=1:5 ; subplot(4,12,i+12) ; topoplot(squeeze(mean(gamma_elecs(:,contrastinds(i),:),1)),eeg.chanlocs,'maplimits',[-.3,.3],'style','both','electrodes','off') ; end ; colormap jet;
for i=1:5 ; subplot(4,12,i+24) ; topoplot(squeeze(mean(alpha_elecs(:,contrastinds(i),:),1)),eeg.chanlocs,'maplimits',[-1.5,1.5],'style','both','electrodes','off') ; end ; colormap jet; 
subplot(4,12,6) ; imagesc([-3,3]) ; h=colorbar; title(h,'dB'); 
subplot(4,12,18) ; imagesc([-.3,.3]) ; h=colorbar; title(h,'\muV (gamma)'); 
subplot(4,12,30) ; imagesc([-1.5,1.5]) ; h=colorbar; title(h,'\muV (alpha/beta)'); 



subplot(1,5,1); 
plot(freqs,squeeze(mean(mean(mersp(:,16,:,times>0 & times<1.75),1),4)),'Color',[0,0,0]); hold on ; 
plot(freqs,squeeze(mean(mean(mersp(:,15,:,times>0 & times<1.75),1),4)),'Color',[0.4,0,0]); 
plot(freqs,squeeze(mean(mean(mersp(:,14,:,times>0 & times<1.75),1),4)),'Color',[0.8,0,0]); 
plot(freqs,squeeze(mean(mean(mersp(:,13,:,times>0 & times<1.75),1),4)),'Color',[1,0.5,0]);
plot(freqs,squeeze(mean(mean(mersp(:,12,:,times>0 & times<1.75),1),4)),'Color',[0.5,0.5,1]); 
hline(0,'k'); 

subplot(1,5,2); 
plot(freqs,squeeze(mean(mean(mersp(:,6,:,times>0 & times<1.75),1),4)),'Color',[0,0,1]); hold on ; 
plot(freqs,squeeze(mean(mean(mersp(:,7,:,times>0 & times<1.75),1),4)),'Color',[0,0.8,.8]); hold on ; 
plot(freqs,squeeze(mean(mean(mersp(:,8,:,times>0 & times<1.75),1),4)),'Color',[0.5,1,0]); hold on ; 
plot(freqs,squeeze(mean(mean(mersp(:,9,:,times>0 & times<1.75),1),4)),'Color',[1,0.5,.1]); hold on ; 
hline(0,'k'); 

subplot(1,5,3); 
plot(freqs,squeeze(mean(mean(mersp(:,2,:,times>0 & times<1.75),1),4)),'Color',[0,0,1]); hold on ; 
plot(freqs,squeeze(mean(mean(mersp(:,3,:,times>0 & times<1.75),1),4)),'Color',[0,0.8,0]); hold on ; 
plot(freqs,squeeze(mean(mean(mersp(:,4,:,times>0 & times<1.75),1),4)),'Color',[0.5,0,0]); hold on ; 
plot(freqs,squeeze(mean(mean(mersp(:,5,:,times>0 & times<1.75),1),4)),'Color',[0,0.5,.1]); hold on ; 
hline(0,'k'); 

subplot(1,5,4); 
plot(freqs,squeeze(mean(mean(mersp(:,1,:,times>0 & times<1.75),1),4)),'Color',[0.8,0,0]); hold on ; 
plot(freqs,squeeze(mean(mean(mersp(:,10,:,times>0 & times<1.75),1),4)),'Color',[0,0.7,0]); hold on ; 
plot(freqs,squeeze(mean(mean(mersp(:,11,:,times>0 & times<1.75),1),4)),'Color',[0,0,0.7]); hold on ; 
hline(0,'k'); 




%{
figure,for i=1:16 ; subplot(4,4,i) ; topoplot(squeeze(mean(gamma_elecs(:,i,:),1)),eeg.chanlocs,'maplimits',[-.3,.3],'style','both','electrodes','off') ; end ; colormap jet;
figure,for i=1:16 ; subplot(4,4,i) ; topoplot(squeeze(mean(alpha_elecs(:,i,:),1)),eeg.chanlocs,'maplimits',[-1.5,1.5],'style','both','electrodes','off') ; end ; colormap jet; 

% make the retinotopic stimuli:
[xg,yg] = meshgrid(-150:150,-150:150) ; 
width = length(xg) ; height = length(yg) ; 
circmask = double(sqrt(xg.^2 + yg.^2) < 350);
fullfield = sqrt(xg.^2 + yg.^2) < 150 ; fullfield = fullfield.*circmask ;
left = (xg<0) .* circmask ; right = (xg>=0) .* circmask ; top = (yg<=0).*circmask ; bot = (yg>0).*circmask ; 
[th,rh] = cart2pol(xg,yg) ; 
wedge1 = (th>0 & th < pi/4).*circmask ; wedge2 = (th >= pi/4 & th < pi/2).*circmask ; 
wedge3 = (th>=pi/2 & th < 3*pi/4).*circmask ; wedge4 = (th>=3*pi/4 & th <= pi).*circmask ; 
wedge5 = (th<=0 & th < -(pi-pi/4)).*circmask ; wedge6 = (th >= -(pi-pi/4) & th < -(pi-pi/2)).*circmask ;  
wedge7 = (th >= -(pi-pi/2) & th < -(pi-3*pi/4)).*circmask ; wedge8 = (th >= -(pi-3*pi/4) & th <= 0).*circmask ; 
allwedge(:,:,1) = wedge1 ; allwedge(:,:,2) = wedge2 ; allwedge(:,:,3) = wedge3 ; allwedge(:,:,4) = wedge4 ; 
allwedge(:,:,5) = wedge5 ; allwedge(:,:,6) = wedge6 ; allwedge(:,:,7) = wedge7 ; allwedge(:,:,8) = wedge8 ; 
ring1 = double(sqrt(xg.^2 + yg.^2) < 400) & double(sqrt(xg.^2 + yg.^2) >=200) ;
ring2 = double(sqrt(xg.^2 + yg.^2) < 200) & double(sqrt(xg.^2 + yg.^2) >=50) ;
ring3 = double(sqrt(xg.^2 + yg.^2) < 50) & double(sqrt(xg.^2 + yg.^2) >=0) ;
allrings(:,:,1) = ring1 ; allrings(:,:,2) = ring2 ; allrings(:,:,3) = ring3 ; 
quad1 = (th>0 & th<pi/2).*circmask ; quad2 = (th>=pi/2 & th<=pi).*circmask ; 
quad3 = (th>=-pi & th<-pi/2).*circmask ; quad4 = (th>=-pi/2 & th<=0).*circmask ; 
allquads(:,:,1) = quad1  ; allquads(:,:,2) = quad2 ; allquads(:,:,3) = quad3 ; allquads(:,:,4) = quad4 ; 
rots = 22.5:22.5:361 ; 
bar1 = (yg<10 & yg>-10).*circmask ; bar2 = imrotate(bar1,rots(1),'nearest','crop').*circmask ; 
bar3 = imrotate(bar1,rots(2),'nearest','crop').*circmask ; bar4 = imrotate(bar1,rots(3),'nearest','crop').*circmask ;
bar5 = imrotate(bar1,rots(4),'nearest','crop').*circmask ; bar6 = imrotate(bar1,rots(5),'nearest','crop').*circmask ;
bar7 = imrotate(bar1,rots(6),'nearest','crop').*circmask ; bar8 = imrotate(bar1,rots(7),'nearest','crop').*circmask ; 
allbars(:,:,1) = bar1 ; allbars(:,:,2) = bar2 ; allbars(:,:,3) = bar3 ; allbars(:,:,4) = bar4 ; 
allbars(:,:,5) = bar5 ; allbars(:,:,6) = bar6 ; allbars(:,:,7) = bar7 ; allbars(:,:,8) = bar8 ; 
clear allstims ; 
allstims(:,:,1) = fullfield ; allstims(:,:,2) = left ; allstims(:,:,3) = right ; allstims(:,:,4) = top ; allstims(:,:,5) = bot ; 
allstims(:,:,6:13) = allwedge ; allstims(:,:,14:17) = allquads ; allstims(:,:,18:20) = allrings ; 
for i=1:8
   allstims(:,:,i+20) = imrotate(allwedge(:,:,i),22.5,'crop') ;  
end
circmask = double(sqrt(xg.^2 + yg.^2) < 70);
clear rhosin
for i=1:200 
    [xg,yg] = meshgrid(-50:50,-50:50) ; 
    [th,rh] = cart2pol(xg,yg) ;
    rhosin = sin(rh+i/16) ;  
end
titles = {'full','left','right','top','bottom','bottom right','bottom left','top left','top right','periphery','fovea'};
retstims = [1,2,3,4,5,14,15,16,17,19,20]; 
figure,for i=1:length(retstims) ; subplot(1,11,i) ; maski = imresize(squeeze(allstims(:,:,retstims(i))).*allstims(:,:,1),[101,101]) ; allmasks(i,:,:) = maski; imagesc(rhosin.*maski); colormap gray; set(gca,'XTickLabel',[],'YTickLabel',[]);  end
allmasks = allmasks > 0.5; 
cd e:/saved; ls 
clipimgs = load('clipimgs.mat') ;clipimgs = clipimgs.clipimgs; 
temp = clipimgs(8,:,:,:) ; clipimgs(8,:,:,:) = clipimgs(9,:,:,:); clipimgs(9,:,:,:) = temp; 

cd C:\Users\butr2901\Desktop\tuning
close all; 
for i=1:11
    %f=figure,
    %subplot(2,3,1); 
    subplottight(5,11,11*1-10+i-1); 
    maski = imresize(squeeze(allstims(:,:,retstims(i))).*allstims(:,:,1),[101,101]) ; imagesc(repmat(mat2gray(rhosin.*maski + 1),[1,1,3])); colormap gray; set(gca,'XTickLabel',[],'YTickLabel',[]);
    %subplot(2,3,4); 
      subplottight(5,11,11*3-10+i-1); 
    imagesc(squeeze(clipimgs(i,:,:,:))); set(gca,'XTickLabel',[],'YTickLabel',[]);
    %subplot(2,3,2);
      subplottight(5,11,11*2-10+i-1); 
    imagesc(times,freqs,squeeze(mean(mersp(:,i,:,:),1)),[-3,3]) ; axis xy ; colormap jet;  set(gca,'XTickLabel',[],'YTickLabel',[]);
    
    %subplot(2,3,5);
      subplottight(5,11,11*4-10+i-1); 
    topoplot(squeeze(mean(gamma_elecs(:,i,:),1)),eeg.chanlocs,'maplimits',[-.3,.3],'style','both','electrodes','off') ;
    %subplot(2,3,6); 
      subplottight(5,11,11*5-10+i-1); 
    topoplot(squeeze(mean(alpha_elecs(:,i,:),1)),eeg.chanlocs,'maplimits',[-1.5,1.5],'style','both','electrodes','off') ;   
   % suptitle(titles{i}); 
    
    %{
    set(f,'Units','normalized');
    set(f,'Position',[0.2 .45 .35 .35]);
    saveas(f,[titles{i},'.svg'],'svg')
    %}
end

for i=1:11 ; subplot(3,4,i) ; imagesc(squeeze(allmasks(i,:,:)>.5)); end

colors = {[1,0,0],[0,1,0],[0,0,1],[0.5,0.5,0],[0,0.5,0.5],[1,0,0],[0,1,0],[0,0,1],[0.5,0.5,0],[1,0,0],[0,1,0]};

titles = {'full','left','right','top','bottom','bottom right','bottom left','top left','top right','periphery','fovea'};
clear allcmasks
figure,for i=1:11 ; subplot(3,4,i) ; maski = squeeze(allmasks(i,:,:)); cmaski(:,:,1) = maski*colors{i}(1); cmaski(:,:,2) = maski*colors{i}(2); cmaski(:,:,3) = maski*colors{i}(3); imagesc(cmaski); allcmasks(i,:,:,:) = cmaski; end 
allcmasks = uint8(allcmasks*255); 
figure,
subplot(2,3,1) ;
imshow(squeeze(allcmasks(1,:,:,:))); title('full');
subplot(2,3,2);
imshow(squeeze(allcmasks(2,:,:,:) + allcmasks(3,:,:,:))); title('left, right'); 
subplot(2,3,3); 
imshow(squeeze(allcmasks(4,:,:,:) + allcmasks(5,:,:,:))); title('top, bottom'); 
subplot(2,3,4); 
imshow(squeeze(allcmasks(6,:,:,:) + allcmasks(7,:,:,:)+allcmasks(8,:,:,:) + allcmasks(9,:,:,:))); title('quadrants'); 
subplot(2,3,5); 
imshow(squeeze(allcmasks(10,:,:,:) + allcmasks(11,:,:,:))); title('periphery, fovea'); 

figure,
subplot(2,6,1) ; 
hz1=freqs>40 & freqs<70; 
hz2 = freqs>8 & freqs<25; 
for i=1:5 
   shadedErrorBar(times,squeeze(mean(mean(mersp(:,i,hz1,:),1),3)),squeeze(std(mean(mersp(:,i,hz1,:),3),0,1))/sqrt(8),{'Color',colors{i}});  hold on ;  
   hline(0,'k') ; vline([0,1.75],'k'); xlim([-.3,2]) ; xlabel('time(s)') ; ylabel('dB'); 
end
subplot(2,6,4) ; for i=1:5 ; plot(1,'Color',colors{i}); hold on; end ; legend(titles(1:5)); 
subplot(2,6,2) ; 
for i=6:9 
   shadedErrorBar(times,squeeze(mean(mean(mersp(:,i,hz1,:),1),3)),squeeze(std(mean(mersp(:,i,hz1,:),3),0,1))/sqrt(8),{'Color',colors{i}});  hold on  ;
   hline(0,'k') ; vline([0,1.75],'k'); xlim([-.3,2]) ;
end
subplot(2,6,5) ; for i=1:4 ; plot(1,'Color',colors{i+5}); hold on; end ; legend(titles(6:9)); 

subplot(2,6,3) ; 
for i=10:11
   shadedErrorBar(times,squeeze(mean(mean(mersp(:,i,hz1,:),1),3)),squeeze(std(mean(mersp(:,i,hz1,:),3),0,1))/sqrt(8),{'Color',colors{i}});  hold on  ;
   hline(0,'k') ; vline([0,1.75],'k'); xlim([-.3,2]) ;
end
subplot(2,6,6) ; for i=1:2 ; plot(1,'Color',colors{i+9}); hold on; end ; legend(titles(10:11)); 

% alpha/beta 
hz = freqs>8 & freqs<25; 
subplot(2,6,7) ; 
for i=1:5 
   shadedErrorBar(times,squeeze(mean(mean(mersp(:,i,hz2,:),1),3)),squeeze(std(mean(mersp(:,i,hz2,:),3),0,1))/sqrt(8),{'Color',colors{i}});  hold on ;  
   hline(0,'k') ; vline([0,1.75],'k'); xlim([-.3,2]) ;
end
subplot(2,6,8) ; 
for i=6:9 
   shadedErrorBar(times,squeeze(mean(mean(mersp(:,i,hz2,:),1),3)),squeeze(std(mean(mersp(:,i,hz2,:),3),0,1))/sqrt(8),{'Color',colors{i}});  hold on  ;
   hline(0,'k') ; vline([0,1.75],'k'); xlim([-.3,2]) ; 
end
subplot(2,6,9) ; 
for i=10:11
   shadedErrorBar(times,squeeze(mean(mean(mersp(:,i,hz2,:),1),3)),squeeze(std(mean(mersp(:,i,hz2,:),3),0,1))/sqrt(8),{'Color',colors{i}});  hold on  ;
   hline(0,'k') ; vline([0,1.75],'k'); xlim([-.3,2]) ; 
end

% get the voxel indices of each different part of the heat map:
titles = {'full','left','right','top','bottom','bottom right','bottom left','top left','top right','periphery','fovea'};
retlabels = {'periphery top left','periphery top right','fovea top left','fovea top right','periphery bottom left','periphery bottom right','fovea bottom left','fovea bottom right'};
retmasks(1,:,:) = squeeze(allmasks(8,:,:).*allmasks(10,:,:)); retmaskinds{1} = find(squeeze(retmasks(1,:,:))~=0); 
retmasks(2,:,:) = squeeze(allmasks(9,:,:).*allmasks(10,:,:)); retmaskinds{2} = find(squeeze(retmasks(2,:,:))~=0); 
retmasks(3,:,:) = squeeze(allmasks(8,:,:).*allmasks(11,:,:)); retmaskinds{3} = find(squeeze(retmasks(3,:,:))~=0); 
retmasks(4,:,:) = squeeze(allmasks(9,:,:).*allmasks(11,:,:)); retmaskinds{4} = find(squeeze(retmasks(4,:,:))~=0); 
retmasks(5,:,:) = squeeze(allmasks(7,:,:).*allmasks(10,:,:)); retmaskinds{5} = find(squeeze(retmasks(5,:,:))~=0); 
retmasks(6,:,:) = squeeze(allmasks(6,:,:).*allmasks(10,:,:)); retmaskinds{6} = find(squeeze(retmasks(6,:,:))~=0); 
retmasks(7,:,:) = squeeze(allmasks(7,:,:).*allmasks(11,:,:)); retmaskinds{7} = find(squeeze(retmasks(7,:,:))~=0); 
retmasks(8,:,:) = squeeze(allmasks(6,:,:).*allmasks(11,:,:)); retmaskinds{8} = find(squeeze(retmasks(8,:,:))~=0); 

% MAKE the heat map out of the stimulus masks (add all together) 
groups = {[2,3],[4,5],[6,7,8,9],[10,11]};
clear gamma_masks alpha_masks 
for i=1:11
    gamma_masks(i,:,:) = squeeze((allmasks(i,:,:)).*squeeze(mean(mean(mean((mersp(:,i,hz1,times>0 & times<1.75)),1),3),4)));
    
    for j=1:8
        sb_gamma_masks(i,j,:,:) = squeeze(allmasks(i,:,:).*squeeze(mean(mean(mersp(j,i,hz1,times>0 & times<1.75),3),4)));
    end
    
    alpha_masks(i,:,:) = squeeze((allmasks(i,:,:)).*squeeze(mean(mean(mean((mersp(:,i,hz2,times>0 & times<1.75)),1),3),4))); 
end



% use z-score as difference from mean for each config
clear groupgamma groupalpha

erodemask = imerode(squeeze(gamma_masks(1,:,:)),strel(ones(3,3,3))); 

for i=1:length(groups) ; groupgamma(i,:,:) = squeeze(mean(gamma_masks(groups{i},:,:),1)); groupalpha(i,:,:) = squeeze(mean(alpha_masks(groups{i},:,:),1));  end
mgroupgamma = squeeze(mean(groupgamma,1)); mgroupgamma = medfilt2(mgroupgamma); zgamma = zscore(unique(mgroupgamma(erodemask~=0))); zgamma(end+1) = -2; zgamma(end+1) = 2; 
%mgroupgamma(mgroupgamma==0) = max(mgroupgamma(:)); % maxgamma = max(mgroupgamma(mgroupgamma~=0)); mingamma = min(mgroupgamma(mgroupgamma ~=0)); 
rgbgamma = ind2rgb(round((mat2gray(mgroupgamma)*255)),parula(255)); 
maskfield = repmat(squeeze(gamma_masks(1,:,:)),[1,1,3]); 
rgbgamma(maskfield==min(maskfield(:))) = 255; 

mgroupalpha = abs(squeeze(mean(groupalpha,1))); mgroupalpha = medfilt2(mgroupalpha); %mgroupalpha(mgroupalpha==0) = max(mgroupalpha(:)); % maxgamma = max(mgroupgamma(mgroupgamma~=0)); mingamma = min(mgroupgamma(mgroupgamma ~=0)); 
rgbalpha = ind2rgb(round((mat2gray(mgroupalpha)*255)),parula(255)); 
rgbalpha(maskfield==min(maskfield(:))) = 255; 


subplot(1,2,1); 
imagesc(rgbgamma); colormap jet; hold on; set(gca,'XTick',[],'YTick',[]);
r= 50;x = 51;y = 51;th = 0:pi/50:2*pi;xunit = r * cos(th) + x;yunit = r * sin(th) + y;
plot(xunit, yunit,'k','LineWidth',4);
r= 17;x = 51;y = 51;th = 0:pi/50:2*pi;xunit = r * cos(th) + x;yunit = r * sin(th) + y;
plot(xunit, yunit,'k','LineWidth',4); 
line([0,101],[51,51],'LineWidth',4,'Color',[0,0,0]);line([51,51],[0,101],'LineWidth',4,'Color',[0,0,0]);

vline(51,'k') ; hline(51,'k'); 

subplot(1,2,2); 
imagesc(rgbalpha) ; colormap jet; hold on; set(gca,'XTick',[],'YTick',[]); 
r= 50;x = 51;y = 51;th = 0:pi/50:2*pi;xunit = r * cos(th) + x;yunit = r * sin(th) + y;
plot(xunit, yunit,'k','LineWidth',4);
r= 17;x = 51;y = 51;th = 0:pi/50:2*pi;xunit = r * cos(th) + x;yunit = r * sin(th) + y;
plot(xunit, yunit,'k','LineWidth',4); 
line([0,101],[51,51],'LineWidth',4,'Color',[0,0,0]);line([51,51],[0,101],'LineWidth',4,'Color',[0,0,0]);

% subject specific tuning fields 
for i=1:length(groups) ; sbfields(:,i,:,:) = squeeze(gamma_masks(groups{i},:,:)); end

cd e:/saved; 
frets = load('all_fmri_rets'); frets = frets.all_fmri_rets; 
fmri_names = {'bottom right','bottom left','top right','top left','left','right','top','bottom','full','fovea','periphery'}; 
titles = {'full','left','right','top','bottom','bottom right','bottom left','top left','top right','periphery','fovea'};
for i=1:length(fmri_names) ; fmri_inds(i) = find(strcmpi(titles{i},fmri_names)) ; end
frets = frets(:,:,:,:,fmri_inds); 

for i=1:11
    fmri_masks(i,:,:) = squeeze((allmasks(i,:,:)).*squeeze(sum(sum(sum(mean(frets(:,:,:,:,i),4) > 0.3)))));    
end
%}


%{

% show the tuning prf using the retinotopic masks 
retnames = {'full','left','right','top','bottom','bottom right','bottom left','top left','top right','periphery','fovea',...
    '1% contrast','5% contrast','15% contrast','35% contrast','100% contrast'};

mtersp = (squeeze(mean(mersp(:,:,:,40:80),4))); 

mquads = squeeze(sum(mtersp(:,6:9,:),2)); 
mlquads = squeeze(sum(mtersp(:,8:9,:),2));
mrquads = squeeze(sum(mtersp(:,6:7,:),2)); 
mhemis = squeeze(sum(mtersp(:,2:3,:),2));
mfields = squeeze(sum(mtersp(:,4:5,:),2)); 
mfovs = squeeze(sum(mtersp(:,10:11,:),2)); 
mfull = squeeze(sum(mtersp(:,1,:),2)); 
subs = 1:8; 

all_configs(1,:,:) = mfull; all_configs(2,:,:) = mhemis; all_configs(3,:,:) = mfields; all_configs(4,:,:) = mfovs; all_configs(5,:,:) = mquads;
labels = {'full','l+r','t+b','f+p','q+q+q+q'};
subplot(2,2,1); 
barwitherr(squeeze(std(mean(all_configs(:,subs,4:12),3),0,2)./sqrt(8)),squeeze(mean(mean(all_configs(:,subs,4:12),2),3))); title('alpha/beta'); 
set(gca,'XTickLabel',labels); 
subplot(2,2,2);
barwitherr(squeeze(std(mean(all_configs(:,subs,15:30),3),0,2)./sqrt(8)),squeeze(mean(mean(all_configs(:,subs,15:30),2),3))); title('gamma');
set(gca,'XTickLabel',labels); 
clear p ts 
for i=1:50
   [~,p(1,i),~,stats] = ttest(mfull(subs,i),mquads(subs,i)); ts(1,i) = stats.tstat;  
   [~,p(2,i),~,stats] = ttest(mfull(subs,i),mhemis(subs,i)); ts(2,i) = stats.tstat;  
   [~,p(3,i),~,stats] = ttest(mfull(subs,i),mfields(subs,i)); ts(3,i) = stats.tstat;  
   [~,p(4,i),~,stats] = ttest(mfull(subs,i),mfovs(subs,i)); ts(4,i) = stats.tstat;  
end
subplot(2,2,3);shadedErrorBar(1:2:100,mean(ts,1),std(ts,0,1)/2); hline(0,'k'); xlabel('frequency(hz)'); ylabel('t-value'); title('mean t-value difference'); vline(30,'r'); 
subplot(2,2,4); imagesc(1:2:100,1:4,p,[0,0.05]); set(gca,'YTick',1:4,'YTickLabel',{'q1+q2+q3+q4','left+right','top+bot','fov+periph'}); xlabel('frequency(hz)'); colorbar; title('p-value for all summed configs'); vline(30,'r');

clear low_h low_p low_ci low_stats low_ts; 
for i=1:4
    [low_h(i),low_p(i),~,low_stats] = ttest(squeeze(mean(all_configs(1,subs,4:12),3)),squeeze(mean(all_configs(i+1,subs,4:12),3))); 
    low_ts(i) = low_stats.tstat; 
    [high_h(i),high_p(i),~,high_stats] = ttest(squeeze(mean(all_configs(1,subs,15:30),3)),squeeze(mean(all_configs(i+1,subs,15:30),3))); 
    high_ts(i) = high_stats.tstat; 
end

all_conts(1,:,:) = mtersp(:,16,:); all_conts(2,:,:) = mtersp(:,15,:); all_conts(3,:,:) = mtersp(:,14,:); all_conts(4,:,:) = mtersp(:,13,:); all_conts(5,:,:) = mtersp(:,12,:);  
labels = {'100%','30%','10%','5%','1%'};
figure,
subplot(2,2,1); 
barwitherr(squeeze(std(mean(all_conts(:,subs,4:12),3),0,2)./sqrt(7)),squeeze(mean(mean(all_conts(:,subs,4:12),2),3))); title('alpha/beta'); 
set(gca,'XTickLabel',labels); 
subplot(2,2,2);
barwitherr(squeeze(std(mean(all_conts(:,subs,15:30),3),0,2)./sqrt(7)),squeeze(mean(mean(all_conts(:,subs,15:30),2),3))); title('gamma');
set(gca,'XTickLabel',labels); 
clear p ts 
for i=1:50
   [~,p(1,i),~,stats] = ttest(squeeze(all_conts(5,:,i)),squeeze(all_conts(4,:,i))); ts(1,i) = stats.tstat;  
   [~,p(2,i),~,stats] = ttest(squeeze(all_conts(5,:,i)),squeeze(all_conts(3,:,i))); ts(2,i) = stats.tstat;  
   [~,p(3,i),~,stats] = ttest(squeeze(all_conts(5,:,i)),squeeze(all_conts(2,:,i))); ts(3,i) = stats.tstat;  
   [~,p(4,i),~,stats] = ttest(squeeze(all_conts(5,:,i)),squeeze(all_conts(1,:,i))); ts(4,i) = stats.tstat;  
end
subplot(2,2,3);shadedErrorBar(1:2:100,mean(ts,1),std(ts,0,1)/2); hline(0,'k'); xlabel('frequency(hz)'); ylabel('t-value'); title('mean t-value difference'); vline(30,'r'); 
subplot(2,2,4); imagesc(1:2:100,1:4,p,[0,0.05]); set(gca,'YTick',1:4,'YTickLabel',{'100% - 30%','100% - 10%','100% - 5%','100% - 1%'}); 
xlabel('frequency(hz)'); colorbar; title('p-value for 100% minus other %contrast'); vline(30,'r');


plot(squeeze(mean(mean(mean(mersp(subs,1:11,5:12,times>0&times<1.5),1),3),4)),squeeze(mean(mean(mean(mersp(subs,1:11,15:35,times>0&times<1.5),1),3),4)),'kd') ; lsline ;
[c,p] = corr(squeeze(mean(mean(mean(mersp(subs,1:11,5:12,times>0&times<1.5),1),3),4))',squeeze(mean(mean(mean(mersp(subs,1:11,15:35,times>0&times<1.5),1),3),4))');
title(['rho=',num2str(c),', p=',num2str(p)]) ; xlabel('alpha/beta'); ylabel('gamma'); 

%}






