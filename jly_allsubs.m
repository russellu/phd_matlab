clear all ; close all
subs = {'russell','camille','ariane','annfred','claudie','gabriel','krisen','Jessily'};
bades = {[16,19,20],[45],[17,22],[],[],[32],[51,52,58],[]};
subcomps ={[5,13,18,19],[3,8],[6,10,14],[3,9],[4,5,9],[5,7,10,15],[3,11,14],[1,14,16]}; 
stims = {'S  1','S  2','S  3','S  4','S  5','S 14','S 15','S 16','S 17','S 19','S 20','S 31','S 32','S 33','S 34','S 35'}; 
cd E:\saved
ret_stersp = load('ret_stersp'); ret_stersp = ret_stersp.ret_stersp; 
nrf_alpha = squeeze(mean(mean(ret_stersp(:,5:9,:),1),2)); 
nrf_gamma = squeeze(mean(mean(ret_stersp(:,15:30,:),1),2)); 
for sb=1:length(subs) 
   
cd(['E:\jly_russ\',subs{sb}]);
%{
vhdrs = dir('*vhdr'); 

for filename = 1:length(vhdrs)
eeg = pop_loadbv('.',vhdrs(filename).name);  % load the raw file
eeg = pop_chanedit(eeg,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
eeg = pop_resample(eeg,250); 
eeg.data = eeg.data - eegfiltfft(eeg.data,eeg.srate,59,61); % new data = old data minus 60Hz
if filename==1; merged = eeg; else; merged = pop_mergeset(merged,eeg); end 
end

%figure,bar(sum(abs(diff(eeg.data,1,2)),2))

merged = pop_interp(merged,bades{sb},'spherical');

pop_saveset(merged,'merged.set'); 
%}
merged = pop_loadset('merged.set'); 
%{
filtdat = eegfiltfft(merged.data,merged.srate,1,125); 
[weights,sphere] = runica(filtdat,'maxsteps',128);
winv = pinv(weights*sphere); 
fullcomps{1} = weights; fullcomps{2} = sphere; save('fullcomps','fullcomps'); 
%figure,for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),merged.chanlocs); title(i); end ; suptitle(subs{sb}); 
%}


fullcomps = load('fullcomps') ; fullcomps = fullcomps.fullcomps; weights = fullcomps{1}; sphere = fullcomps{2}; 
comp_inds = load('comp_inds'); comp_inds = comp_inds.comp_inds; 
winv = pinv(weights*sphere); 
bads = zeros(1,64) ; bads(comp_inds(1:5)) = 1; bads = find(bads==0); 
icadat = weights*sphere*merged.data; 
icadat(bads,:) = 0;
denoised_channels = winv*icadat; 

comp_merged = merged; comp_merged.data = weights*sphere*merged.data; 


gamma_denoised = (eegfiltfft(denoised_channels,merged.srate,20,120)); 
gamma_merged = merged; gamma_merged.data = gamma_denoised;
clear merspsorts 
clear mstersp; 
for st=1:length(stims)
    allep = pop_epoch(comp_merged,{stims{st}},[-.5,2.5]);
    noiseep = pop_epoch(gamma_merged,{stims{st}},[-.5,1.5]);
    [sv,si] = sort(mean(std(noiseep.data(comp_inds(1:5),:,:),0,2),1),'descend'); 
    merspsorts{st} = si; 
    for i=1:5
        [mstersp(st,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(comp_inds(i),:,si(8:end))),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',80,'baseline',0,'verbose','off','timesout',100) ; 
    end
end
%figure,for i=1:16 ; subplot(4,4,i) ; imagesc(squeeze(mean(mstersp(i,:,:,:),2)),[-3,3]); axis xy ; colormap jet; end ; suptitle(subs{sb}); 
allersp(sb,:,:,:) = squeeze(mean(mstersp,2)); 


alphabeta_denoised = abs(eegfiltfft(denoised_channels,merged.srate,10,16)); 
gamma_denoised = abs(eegfiltfft(denoised_channels,merged.srate,35,60)); 
gamma_merged = merged; gamma_merged.data = gamma_denoised;
alphabeta_merged = merged; alphabeta_merged.data = alphabeta_denoised;

for st=1:length(stims)
    ep_gamma = pop_epoch(gamma_merged,{stims{st}},[-.5,2.5]);
    ep_alphabeta = pop_epoch(alphabeta_merged,{stims{st}},[-.5,2.5]);
    
    [sv,si] = sort(squeeze(mean(std(ep_gamma.data,0,2),1)),'descend'); 
    
    mgamma(sb,st,:) = squeeze(mean(mean(ep_gamma.data(:,ep_gamma.times>0 & ep_gamma.times<2000,merspsorts{st}(8:end)),2)-mean(ep_gamma.data(:,ep_gamma.times<0,merspsorts{st}(8:end)),2),3)); 
    malphabeta(sb,st,:) = squeeze(mean(mean(ep_alphabeta.data(:,ep_alphabeta.times>0 & ep_alphabeta.times<2000,merspsorts{st}(8:end)),2)-mean(ep_alphabeta.data(:,ep_alphabeta.times<0,merspsorts{st}(8:end)),2),3)); 
end


%{
fullcomps = load('fullcomps') ; fullcomps = fullcomps.fullcomps; weights = fullcomps{1}; sphere = fullcomps{2}; 
ica_merged = merged; ica_merged.data = weights*sphere*merged.data; 

allep = pop_epoch(ica_merged,stims,[-.5,2.25]);

for i=1:64; disp(i); 
    for j=1:size(allep.data,3)
        [mstersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,j)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
            'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',80,'baseline',0,'verbose','off','timesout',100) ; 
    end
end

mgamma = squeeze(mean(mstersp(:,:,15:30,:),3)); 
malpha = squeeze(mean(mstersp(:,:,5:9,:),3)); 

for i=1:size(mgamma,1)
    for j=1:size(mgamma,2)
        gammacorrs(i,j) = corr2(nrf_gamma,squeeze(mgamma(i,j,:))); 
        alphacorrs(i,j) = corr2(nrf_alpha,squeeze(malpha(i,j,:))); 
    end
end


[sv,si] = sort(mean(gammacorrs,2) + mean(alphacorrs,2),'descend'); 
figure, for i=1:length(si) ; subplot(5,13,i) ; imagesc(squeeze(mean(mstersp(si(i),:,:,:),2)),[-6,2]) ; axis xy ; colormap jet; end ; suptitle(subs{sb}); 

comp_inds = si; save('comp_inds','comp_inds'); 
%}
%{
for i=1:64
    [mstersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
        'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',80,'baseline',0,'verbose','off','timesout',100) ; 
end
allstersp(sb,:,:) = squeeze(mean(mstersp(subcomps{sb},:,:),1));
%}
%figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mstersp(i,:,:)),[-3,3]) ; axis xy ; colormap jet; title(i); end; suptitle(subs{sb}); 

end


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
circmask = double(sqrt(xg.^2 + yg.^2) < 50);
clear rhosin
for i=1:200 
    [xg,yg] = meshgrid(-50:50,-50:50) ; 
    [th,rh] = cart2pol(xg,yg) ;
    rhosin = sin(rh+i/16) ;  
end
titles = {'full','left','right','top','bottom','bottom right','bottom left','top left','top right','periphery','fovea'};
retstims = [1,2,3,4,5,14,15,16,17,19,20]; 
clear allmasks 
figure,for i=1:length(retstims) ; subplot(1,11,i) ; maski = imresize(squeeze(allstims(:,:,retstims(i))).*allstims(:,:,1),[101,101]) ; allmasks(i,:,:) = maski; imagesc(rhosin.*maski); colormap gray; set(gca,'XTickLabel',[],'YTickLabel',[]);  end
allmasks = allmasks >=.5; 

retnames = {'full','left','right','top','bottom','bottom right','bottom left','top left','top right','periphery','fovea','1%','5%','10%','30%','100%'};

retino_colors = [[127,65,69]/255;[189,61,58]/255;[63,105,170]/255;[213,174,65]/255;[118,111,87]/255;[228,122,46]/255;[190,158,201]/255;[241,234,127]/255;[5,110,109]/255;[72,81,103]/255;[209,184,148]/255]; 
repmasks = double(repmat(allmasks,[1,1,1,3])); 
for i=1:11 
    repmasks(i,:,:,1) = double(repmasks(i,:,:,1)).*retino_colors(i,1); 
    repmasks(i,:,:,2) = double(repmasks(i,:,:,2)).*retino_colors(i,2); 
    repmasks(i,:,:,3) = double(repmasks(i,:,:,3)).*retino_colors(i,3); 
end

stim_masks(1,:,:,:) = repmasks(1,:,:,:); 
stim_masks(2,:,:,:) = repmasks(2,:,:,:) + repmasks(3,:,:,:); 
stim_masks(3,:,:,:) = repmasks(4,:,:,:) + repmasks(5,:,:,:); 
stim_masks(4,:,:,:) = repmasks(6,:,:,:) + repmasks(7,:,:,:) + repmasks(8,:,:,:) + repmasks(9,:,:,:); 
stim_masks(5,:,:,:) = repmasks(10,:,:,:) + repmasks(11,:,:,:);
stim_masks(stim_masks==0) = 255; 
config_titles = {'full','left, right','upper, lower','quadrants','fovea, periphery'};

% retinotopic correlation
subplot(2,2,3); 
plot(squeeze(mean(mean(mgamma(:,ret_inds,all_elecs),1),3)),squeeze(mean(mean(malphabeta(:,ret_inds,all_elecs),1),3)),'o'); lsline; 
[rho, p] = corr(squeeze(mean(mean(mgamma(:,ret_inds,all_elecs),1),3))',squeeze(mean(mean(malphabeta(:,ret_inds,all_elecs),1),3))'); 
title([format_rho(rho),', ',format_p(p)]); 

subplot(2,2,4); 
% contrast correlation 
plot(squeeze(mean(mean(mgamma(:,contrast_inds,all_elecs),1),3)),squeeze(mean(mean(malphabeta(:,contrast_inds,all_elecs),1),3)),'o'); lsline; 
[rho, p] = corr(squeeze(mean(mean(mgamma(:,contrast_inds,all_elecs),1),3))',squeeze(mean(mean(malphabeta(:,contrast_inds,all_elecs),1),3))'); 
title([format_rho(rho),', ',format_p(p)]); 

subplot(4,4,1); 
imshow(squeeze(stim_masks(1,:,:,:))) ; set(gca,'XTick',[],'YTick',[]) ; title(config_titles{1});
subplot(4,4,2); 
imshow(squeeze(stim_masks(2,:,:,:))) ; set(gca,'XTick',[],'YTick',[]) ; title(config_titles{2});
subplot(4,4,3); 
imshow(squeeze(stim_masks(4,:,:,:))) ; set(gca,'XTick',[],'YTick',[]) ; title(config_titles{4});


stim_inds = [1,2,3,7,8,6,9];

figure
for i=1:length(stim_inds) 
    subplot(6,12,i) ;
    maski = imresize(squeeze(allstims(:,:,retstims(stim_inds(i)))).*allstims(:,:,1),[101,101]) ; 
    allmasks(i,:,:) = maski; imagesc(rhosin.*maski); 
    colormap gray; set(gca,'XTickLabel',[],'YTickLabel',[]);  
    title(titles{stim_inds(i)}); 
end
figure,
for i=1:length(stim_inds)
    subplot(6,12,i);
    imagesc(times,freqs,squeeze(mean(allersp(1:6,stim_inds(i),:,:),1)),[-3,3]); axis xy ; colormap jet; 
    if i==1 ; xlabel('time(s)'); ylabel('frequency(hz)'); end
end
subplot(6,12,12) ; imagesc([-3,3]); h = colorbar ; title(h,'dB'); 
for i=1:length(stim_inds)
   subplot(6,12,i+12);
   topoplot(squeeze(mean(mgamma(1:6,stim_inds(i),:),1)),merged.chanlocs,'maplimits',[-.25,.25]);
end
subplot(6,12,24) ; imagesc([-.25,.25]); h = colorbar ; title(h,'\muV'); 
for i=1:length(stim_inds)
   subplot(6,12,i+24);
   topoplot(squeeze(mean(malphabeta(1:6,stim_inds(i),:),1)),merged.chanlocs,'maplimits',[-1,1]);
end
subplot(6,12,36) ; imagesc([-1,1]); h = colorbar ; title(h,'\muV'); 



right_elecs = [20,54,58,26,59,63,64,31];
left_elecs = [52,19,56,24,57,60,61,29];
all_elecs = [right_elecs,left_elecs,25,63,30];

subplot(1,2,1);topoplot(zeros(1,64),merged.chanlocs,'emarker2',{right_elecs,'o','r'}); hold on; 
topoplot(zeros(1,64),merged.chanlocs,'emarker2',{left_elecs,'o','b'}); 


figure,



subplot(4,4,5) ; 
b=bar(1,squeeze(mean(mean(mgamma(1:6,stim_inds(2),left_elecs),1),3))); hold on ; 
b.FaceColor = retino_colors(stim_inds(2),:); 
errorbar(1,squeeze(mean(mean(mgamma(1:6,stim_inds(2),left_elecs),1),3)),squeeze(mean(std(mgamma(1:6,stim_inds(2),left_elecs),0,1),3))/sqrt(6),'k.'); 
b=bar(2,squeeze(mean(mean(mgamma(1:6,stim_inds(3),left_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(3),:); 
errorbar(2,squeeze(mean(mean(mgamma(1:6,stim_inds(3),left_elecs),1),3)),squeeze(mean(std(mgamma(1:6,stim_inds(3),left_elecs),0,1),3))/sqrt(6),'k.'); 
[h,p,ci,stats] = ttest(squeeze(mean(mgamma(1:6,stim_inds(2),left_elecs),3)),squeeze(mean(mgamma(1:6,stim_inds(3),left_elecs),3)));
text(1.5,0.25,format_p(p)); 

b=bar(4,squeeze(mean(mean(mgamma(1:6,stim_inds(4),left_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(4),:); 
errorbar(4,squeeze(mean(mean(mgamma(1:6,stim_inds(4),left_elecs),1),3)),squeeze(mean(std(mgamma(1:6,stim_inds(4),left_elecs),0,1),3))/sqrt(6),'k.'); 
b=bar(5,squeeze(mean(mean(mgamma(1:6,stim_inds(6),left_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(6),:); 
errorbar(5,squeeze(mean(mean(mgamma(1:6,stim_inds(6),left_elecs),1),3)),squeeze(mean(std(mgamma(1:6,stim_inds(6),left_elecs),0,1),3))/sqrt(6),'k.'); 
[h,p,ci,stats] = ttest(squeeze(mean(mgamma(1:6,stim_inds(4),left_elecs),3)),squeeze(mean(mgamma(1:6,stim_inds(6),left_elecs),3)));
text(4.5,0.25,format_p(p)); 

b=bar(7,squeeze(mean(mean(mgamma(1:6,stim_inds(5),left_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(5),:); 
errorbar(7,squeeze(mean(mean(mgamma(1:6,stim_inds(5),left_elecs),1),3)),squeeze(mean(std(mgamma(1:6,stim_inds(5),left_elecs),0,1),3))/sqrt(6),'k.'); 
b=bar(8,squeeze(mean(mean(mgamma(1:6,stim_inds(7),left_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(7),:); 
errorbar(8,squeeze(mean(mean(mgamma(1:6,stim_inds(7),left_elecs),1),3)),squeeze(mean(std(mgamma(1:6,stim_inds(7),left_elecs),0,1),3))/sqrt(6),'k.'); 
[h,p,ci,stats] = ttest(squeeze(mean(mgamma(1:6,stim_inds(5),left_elecs),3)),squeeze(mean(mgamma(1:6,stim_inds(7),left_elecs),3)));
text(7.5,0.15,format_p(p)); 
set(gca,'XTick',[1.5,4.5,7.5],'XTickLabel',{'L vs R','BL vs BR','TL vs TR'}); ylim([0,0.3]); title('left electrodes (gamma)'); 


subplot(4,4,6) ; 
b=bar(1,squeeze(mean(mean(mgamma(1:6,stim_inds(2),right_elecs),1),3))); hold on ; 
b.FaceColor = retino_colors(stim_inds(2),:); 
errorbar(1,squeeze(mean(mean(mgamma(1:6,stim_inds(2),right_elecs),1),3)),squeeze(mean(std(mgamma(1:6,stim_inds(2),right_elecs),0,1),3))/sqrt(6),'k.'); 
b=bar(2,squeeze(mean(mean(mgamma(1:6,stim_inds(3),right_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(3),:); 
errorbar(2,squeeze(mean(mean(mgamma(1:6,stim_inds(3),right_elecs),1),3)),squeeze(mean(std(mgamma(1:6,stim_inds(3),right_elecs),0,1),3))/sqrt(6),'k.'); 
[h,p,ci,stats] = ttest(squeeze(mean(mgamma(1:6,stim_inds(2),right_elecs),3)),squeeze(mean(mgamma(1:6,stim_inds(3),right_elecs),3)));
text(1.5,0.25,format_p(p)); 

b=bar(4,squeeze(mean(mean(mgamma(1:6,stim_inds(4),right_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(4),:); 
errorbar(4,squeeze(mean(mean(mgamma(1:6,stim_inds(4),right_elecs),1),3)),squeeze(mean(std(mgamma(1:6,stim_inds(4),right_elecs),0,1),3))/sqrt(6),'k.'); 
b=bar(5,squeeze(mean(mean(mgamma(1:6,stim_inds(6),right_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(6),:); 
errorbar(5,squeeze(mean(mean(mgamma(1:6,stim_inds(6),right_elecs),1),3)),squeeze(mean(std(mgamma(1:6,stim_inds(6),right_elecs),0,1),3))/sqrt(6),'k.'); 
[h,p,ci,stats] = ttest(squeeze(mean(mgamma(1:6,stim_inds(4),right_elecs),3)),squeeze(mean(mgamma(1:6,stim_inds(6),right_elecs),3)));
text(4.5,0.25,format_p(p)); 

b=bar(7,squeeze(mean(mean(mgamma(1:6,stim_inds(5),right_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(5),:); 
errorbar(7,squeeze(mean(mean(mgamma(1:6,stim_inds(5),right_elecs),1),3)),squeeze(mean(std(mgamma(1:6,stim_inds(5),right_elecs),0,1),3))/sqrt(6),'k.'); 
b=bar(8,squeeze(mean(mean(mgamma(1:6,stim_inds(7),right_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(7),:); 
errorbar(8,squeeze(mean(mean(mgamma(1:6,stim_inds(7),right_elecs),1),3)),squeeze(mean(std(mgamma(1:6,stim_inds(7),right_elecs),0,1),3))/sqrt(6),'k.'); 
[h,p,ci,stats] = ttest(squeeze(mean(mgamma(1:6,stim_inds(5),right_elecs),3)),squeeze(mean(mgamma(1:6,stim_inds(7),right_elecs),3)));
text(7.5,0.15,format_p(p)); 
set(gca,'XTick',[1.5,4.5,7.5],'XTickLabel',{'L vs R','BL vs BR','TL vs TR'}); ylim([0,0.3]); title('right electrodes (gamma)'); 



subplot(4,4,7) ; 
b=bar(1,squeeze(mean(mean(malphabeta(1:6,stim_inds(2),left_elecs),1),3))); hold on ; 
b.FaceColor = retino_colors(stim_inds(2),:); 
errorbar(1,squeeze(mean(mean(malphabeta(1:6,stim_inds(2),left_elecs),1),3)),squeeze(mean(std(malphabeta(1:6,stim_inds(2),left_elecs),0,1),3))/sqrt(6),'k.'); 
b=bar(2,squeeze(mean(mean(malphabeta(1:6,stim_inds(3),left_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(3),:); 
errorbar(2,squeeze(mean(mean(malphabeta(1:6,stim_inds(3),left_elecs),1),3)),squeeze(mean(std(malphabeta(1:6,stim_inds(3),left_elecs),0,1),3))/sqrt(6),'k.'); 
[h,p,ci,stats] = ttest(squeeze(mean(malphabeta(1:6,stim_inds(2),left_elecs),3)),squeeze(mean(malphabeta(1:6,stim_inds(3),left_elecs),3)));
text(1.5,0.1,format_p(p)); 

b=bar(4,squeeze(mean(mean(malphabeta(1:6,stim_inds(4),left_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(4),:); 
errorbar(4,squeeze(mean(mean(malphabeta(1:6,stim_inds(4),left_elecs),1),3)),squeeze(mean(std(malphabeta(1:6,stim_inds(4),left_elecs),0,1),3))/sqrt(6),'k.'); 
b=bar(5,squeeze(mean(mean(malphabeta(1:6,stim_inds(6),left_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(6),:); 
errorbar(5,squeeze(mean(mean(malphabeta(1:6,stim_inds(6),left_elecs),1),3)),squeeze(mean(std(malphabeta(1:6,stim_inds(6),left_elecs),0,1),3))/sqrt(6),'k.'); 
[h,p,ci,stats] = ttest(squeeze(mean(malphabeta(1:6,stim_inds(4),left_elecs),3)),squeeze(mean(malphabeta(1:6,stim_inds(6),left_elecs),3)));
text(4.5,0.1,format_p(p)); 

b=bar(7,squeeze(mean(mean(malphabeta(1:6,stim_inds(5),left_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(5),:); 
errorbar(7,squeeze(mean(mean(malphabeta(1:6,stim_inds(5),left_elecs),1),3)),squeeze(mean(std(malphabeta(1:6,stim_inds(5),left_elecs),0,1),3))/sqrt(6),'k.'); 
b=bar(8,squeeze(mean(mean(malphabeta(1:6,stim_inds(7),left_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(7),:); 
errorbar(8,squeeze(mean(mean(malphabeta(1:6,stim_inds(7),left_elecs),1),3)),squeeze(mean(std(malphabeta(1:6,stim_inds(7),left_elecs),0,1),3))/sqrt(6),'k.'); 
[h,p,ci,stats] = ttest(squeeze(mean(malphabeta(1:6,stim_inds(5),left_elecs),3)),squeeze(mean(malphabeta(1:6,stim_inds(7),left_elecs),3)));
text(7.5,0.1,format_p(p)); 
set(gca,'XTick',[1.5,4.5,7.5],'XTickLabel',{'L vs R','BL vs BR','TL vs TR'}); title('left electrodes (alpha/beta)'); 


subplot(4,4,8) ; 
b=bar(1,squeeze(mean(mean(malphabeta(1:6,stim_inds(2),right_elecs),1),3))); hold on ; 
b.FaceColor = retino_colors(stim_inds(2),:); 
errorbar(1,squeeze(mean(mean(malphabeta(1:6,stim_inds(2),right_elecs),1),3)),squeeze(mean(std(malphabeta(1:6,stim_inds(2),right_elecs),0,1),3))/sqrt(6),'k.'); 
b=bar(2,squeeze(mean(mean(malphabeta(1:6,stim_inds(3),right_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(3),:); 
errorbar(2,squeeze(mean(mean(malphabeta(1:6,stim_inds(3),right_elecs),1),3)),squeeze(mean(std(malphabeta(1:6,stim_inds(3),right_elecs),0,1),3))/sqrt(6),'k.'); 
[h,p,ci,stats] = ttest(squeeze(mean(malphabeta(1:6,stim_inds(2),right_elecs),3)),squeeze(mean(malphabeta(1:6,stim_inds(3),right_elecs),3)));
text(1.5,0.1,format_p(p)); 

b=bar(4,squeeze(mean(mean(malphabeta(1:6,stim_inds(4),right_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(4),:); 
errorbar(4,squeeze(mean(mean(malphabeta(1:6,stim_inds(4),right_elecs),1),3)),squeeze(mean(std(malphabeta(1:6,stim_inds(4),right_elecs),0,1),3))/sqrt(6),'k.'); 
b=bar(5,squeeze(mean(mean(malphabeta(1:6,stim_inds(6),right_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(6),:); 
errorbar(5,squeeze(mean(mean(malphabeta(1:6,stim_inds(6),right_elecs),1),3)),squeeze(mean(std(malphabeta(1:6,stim_inds(6),right_elecs),0,1),3))/sqrt(6),'k.'); 
[h,p,ci,stats] = ttest(squeeze(mean(malphabeta(1:6,stim_inds(4),right_elecs),3)),squeeze(mean(malphabeta(1:6,stim_inds(6),right_elecs),3)));
text(4.5,0.1,format_p(p)); 

b=bar(7,squeeze(mean(mean(malphabeta(1:6,stim_inds(5),right_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(5),:); 
errorbar(7,squeeze(mean(mean(malphabeta(1:6,stim_inds(5),right_elecs),1),3)),squeeze(mean(std(malphabeta(1:6,stim_inds(5),right_elecs),0,1),3))/sqrt(6),'k.'); 
b=bar(8,squeeze(mean(mean(malphabeta(1:6,stim_inds(7),right_elecs),1),3))); 
b.FaceColor = retino_colors(stim_inds(7),:); 
errorbar(8,squeeze(mean(mean(malphabeta(1:6,stim_inds(7),right_elecs),1),3)),squeeze(mean(std(malphabeta(1:6,stim_inds(7),right_elecs),0,1),3))/sqrt(6),'k.'); 
[h,p,ci,stats] = ttest(squeeze(mean(malphabeta(1:6,stim_inds(5),right_elecs),3)),squeeze(mean(malphabeta(1:6,stim_inds(7),right_elecs),3)));
text(7.5,0.1,format_p(p)); 
set(gca,'XTick',[1.5,4.5,7.5],'XTickLabel',{'L vs R','BL vs BR','TL vs TR'}); title('right electrodes (alpha/beta)'); 



%{
figure,imagesc(squeeze(mean(allstersp,1)),[-2,2]) ; axis xy ; colormap jet;
cd E:\saved
ret_stersp = allstersp; save('ret_stersp','ret_stersp'); 
%}
%{
figure,
stim_inds = [1,2,3,4,5,10,11];
for i=1:length(stim_inds) ; subplot(2,7,i) ; topoplot(squeeze(mean(malphabeta(1:6,stim_inds(i),:),1)),merged.chanlocs,'maplimits',[-2,2]); title(retnames{stim_inds(i)}); end
for i=1:length(stim_inds) ; subplot(2,7,i+7) ; topoplot(squeeze(mean(mgamma(1:6,stim_inds(i),:),1)),merged.chanlocs,'maplimits',[-.3,.3]); title(retnames{stim_inds(i)}); end
suptitle('retinotopy'); 
figure,
cont_inds = 16:-1:12; 
for i=1:length(cont_inds) ; subplot(2,5,i) ; topoplot(squeeze(mean(malphabeta(1:6,cont_inds(i),:),1)),merged.chanlocs,'maplimits',[-1.5,1.5]); title(retnames{cont_inds(i)}); end
for i=1:length(cont_inds) ; subplot(2,5,i+5) ; topoplot(squeeze(mean(mgamma(1:6,cont_inds(i),:),1)),merged.chanlocs,'maplimits',[-.3,.3]); title(retnames{cont_inds(i)}); end
suptitle('contrast'); 

topoplot(zeros(1,64),merged.chanlocs,'electrodes','numbers'); 
right_elecs = [20,54,58,26,59,63,64,31];
left_elecs = [52,19,56,24,57,60,61,29];
all_elecs = [right_elecs,left_elecs,25,63,30];

mgamma_right = squeeze(mean(mgamma(1:6,:,right_elecs),3)); 
mgamma_left = squeeze(mean(mgamma(1:6,:,left_elecs),3)); 
malphabeta_right = squeeze(mean(malphabeta(1:6,:,right_elecs),3)); 
malphabeta_left = squeeze(mean(malphabeta(1:6,:,left_elecs),3)); 

for i=1:length(stim_inds) 
    subplot(3,7,i); 
    b = bar(1,squeeze(mean(mgamma_left(:,stim_inds(i)),1))); hold on; 
    b.FaceColor = [0,0,1]; 
    errorbar(1,squeeze(mean(mgamma_left(:,stim_inds(i)),1)),squeeze(std(mgamma_left(:,stim_inds(i)),0,1))/sqrt(6),'k.');

    b = bar(2,squeeze(mean(mgamma_right(:,stim_inds(i)),1))); 
    b.FaceColor = [1,0,1]; 
    errorbar(2,squeeze(mean(mgamma_right(:,stim_inds(i)),1)),squeeze(std(mgamma_right(:,stim_inds(i)),0,1))/sqrt(6),'k.');

    [h,p,ci,stats] = ttest(mgamma_left(:,stim_inds(i)),mgamma_right(:,stim_inds(i))); 
    set(gca,'XTick',[1,2],'XTickLabel',{'left','right'}); xlabel('electrode ROI'); 
    title(format_p(p)); 

end

for i=1:length(stim_inds) 
    subplot(3,7,i+7); 
    b = bar(1,squeeze(mean(malphabeta_left(:,stim_inds(i)),1))); hold on; 
    b.FaceColor = [0,0,1]; 
    errorbar(1,squeeze(mean(malphabeta_left(:,stim_inds(i)),1)),squeeze(std(malphabeta_left(:,stim_inds(i)),0,1))/sqrt(6),'k.');

    b = bar(2,squeeze(mean(malphabeta_right(:,stim_inds(i)),1))); 
    b.FaceColor = [1,0,1]; 
    errorbar(2,squeeze(mean(malphabeta_right(:,stim_inds(i)),1)),squeeze(std(malphabeta_right(:,stim_inds(i)),0,1))/sqrt(6),'k.');
    set(gca,'XTick',[1,2],'XTickLabel',{'left','right'}); %xlabel('electrode ROI hemi'); 
    [h,p,ci,stats] = ttest(malphabeta_left(:,stim_inds(i)),malphabeta_right(:,stim_inds(i))); 
    title(format_p(p)); 

end

CONTRAST_COLORS = {[194,82,60]/255,[230,142,28]/255,[123,237,0]/255,[30,144,148]/255,[11,44,122]/255}; 
subplot(3,7,15); 
for i=1:length(cont_inds)
    b = bar(i,squeeze(mean(mean(mgamma(:,cont_inds(i),all_elecs),1),3))); hold on ;
    b.FaceColor = CONTRAST_COLORS{i}; 
    errorbar(i,squeeze(mean(mean(mgamma(:,cont_inds(i),all_elecs),1),3)),squeeze(std(mean(mgamma(:,cont_inds(i),all_elecs),3),0,1))/sqrt(6),'k.'); 
end

subplot(3,7,16); 
for i=1:length(cont_inds)
    b = bar(i,squeeze(mean(mean(malphabeta(:,cont_inds(i),all_elecs),1),3))); hold on ;
    b.FaceColor = CONTRAST_COLORS{i}; 
    errorbar(i,squeeze(mean(mean(malphabeta(:,cont_inds(i),all_elecs),1),3)),squeeze(std(mean(malphabeta(:,cont_inds(i),all_elecs),3),0,1))/sqrt(6),'k.'); 
end


%}


















