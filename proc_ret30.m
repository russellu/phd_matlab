clear all ; close all ; 
cd C:\shared\orientation
discr = dir('russ_rotation*vhdr') ; 
for i=1:length(discr)
   EEG = pop_loadbv('.',discr(i).name) ; 
   EEG = pop_resample(EEG,256) ; 
   if i==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG) ; end
end
merged = pop_chanedit(merged,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
temp = merged.data(1:32,:) ; merged.data(1:32,:) = merged.data(33:64,:) ; merged.data(33:64,:) = temp ; 
mergefilt = merged ; mergefilt.data = eegfiltfft(merged.data,merged.srate,1,70) ; 
ica = pop_runica(mergefilt,'runica') ; ica.icaact = ica.icaweights*ica.icasphere*merged.data ; 
ica.data = ica.icaact ; 
stims = {'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8','S  9','S 10','S 11','S 12','S 13','S 14','S 15','S 16','S 17','S 18','S 19','S 20','S 21','S 22','S 23','S 24'} ;
clear allersp ; 
for s=1%:length(stims) ; disp(s) ; 
    allep = pop_epoch(ica,{stims{s}},[-2,31]) ;
    for i=1:64 ;
            [allersp(s,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',256,'baseline',0,'verbose','off') ; 
                
    end
end
for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(allersp(:,i,:,:),1)),[-3,3]) ; title(i) ; end

goodcs = [26,28,49] ; 
clear tersp
for s=1%:length(stims) ; disp(s) ; 
    allep = pop_epoch(ica,{stims{s}},[-2,32]) ;
    for i=1:length(goodcs) ;
        for j=1:46
            [tersp(s,i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(goodcs(i),:,j)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',NaN,'verbose','off') ; 
        end       
    end
end
tersp = squeeze(tersp - repmat(mean(tersp(:,:,:,:,times>30 | times<0),5),[1,1,1,1,200])) ; 
for i=1:3 ; for j=1:46 ; medtersp(i,j,:,:) = medfilt2(squeeze(tersp(i,j,:,:)),[3,5]) ; end ; end 
figure,for i=1:3 ; subplot(2,2,i) ; imagesc(medfilt2(squeeze(mean(tersp(i,:,:,:),2))),[-3,3]) ; end
mtersp = squeeze(mean(mean(tersp))) ; 
imagesc(medfilt2(mtersp,[3,5]),[-2,2]) ;
figure,
for i=1:3 ; subplot(2,2,i) ; 
     shadedErrorBar([],squeeze(mean(mean(medtersp(i,:,freqs>40 & freqs<65,times>2 & times<28),2),3)),...
             squeeze(std(mean(medtersp(i,:,freqs>40 & freqs<65,times>2 & times<28),3),0,2))./sqrt(46),{'r'}) ; 
         hold on ; 
     shadedErrorBar([],squeeze(mean(mean(medtersp(i,:,freqs>10 & freqs<16,times>2 & times<28),2),3)),...
     squeeze(std(mean(medtersp(i,:,freqs>10 & freqs<16,times>2 & times<28),3),0,2))./sqrt(46),{'b'}) ; hline(0,'k') ; 
end

plot(squeeze(mean(mean(medtersp(2,:,freqs>45&freqs<60,times>2&times<28),2),3)),...
     squeeze(mean(mean(medtersp(2,:,freqs>8&freqs<18,times>2&times<28),2),3)),'o') ; lsline 


bads = zeros(1,64) ; bads(goodcs) = 1 ; bads = find(bads==0) ; 
acts = ica.icaact ; acts(bads,:) = 0 ; 
invacts = pinv(ica.icaweights*ica.icasphere)*acts ; 
ica.data = invacts ; postelecs = [60,29,30,31,64,63,62,61,23,56,24,57,25,58,26,59,27,28,32] ; 

epdata = pop_epoch(ica,{'S  1'},[-2,32]) ; 
freqs = 1:4:100 ; allfreqs = zeros(size(epdata.data,1),length(freqs),size(epdata.data,2),size(epdata.data,3)) ; 
for f=1:length(freqs)
    for j=1:64
        filtj = eegfiltfft(squeeze(epdata.data(j,:,:))',ica.srate,freqs(f)-3,freqs(f)+3) ;
        allfreqs(j,f,:,:) = filtj' ; 
    end
end
allfreqs = abs(allfreqs) ; 
for i=1:64
    for j=1:25
        dataij = squeeze(allfreqs(i,j,:,:)) ; 
        allfreqs(i,j,:,:) = dataij - repmat(mean(dataij(epdata.times<0 | epdata.times>30000,:),1),[8704,1]) ; 
    end
end

mfreqs = squeeze(mean(allfreqs,4)) ; 
smthfreqs = zeros(size(mfreqs)) ; for i=1:64 ; smthfreqs(i,:,:) = imfilter(squeeze(mfreqs(i,:,:)),fspecial('gaussian',[1,125],125)) ; end
icount = 1 ; clear xg ; 
for i=1:10:8700; disp(i) ; 
    %subplot(8,10,icount) ; 
    [~,xg(:,:,icount)] = topoplot(squeeze(mean(mean(smthfreqs(:,freqs>10 & freqs<20,i:i+10),2),3)),ica.chanlocs,'maplimits',[-.05,.05],'noplot','on') ; 
    [f,xgg(:,:,icount)] = topoplot(squeeze(mean(mean(smthfreqs(:,freqs>40 & freqs<65,i:i+10),2),3)),ica.chanlocs,'maplimits',[-.05,.05],'noplot','on') ; 
    icount = icount + 1 ; 
end
xg2 = xg ; for i=1:length(xg) ; xg2(:,:,i) = flipud(xg(:,:,i)) ; end ; xgg2 = xgg ; for i=1:length(xgg) ; xgg2(:,:,i) = flipud(xgg(:,:,i)) ; end 
labs = {ica.chanlocs.labels} ; 
plot(squeeze(mean(smthfreqs(:,freqs>40 & freqs<65,:),2))') ; 
subplot(1,2,1) ; 
imagesc(squeeze(mean(smthfreqs(:,freqs> 60& freqs<80,:),2))) ; set(gca,'YTick',1:64,'YTickLabel',labs) ;
subplot(1,2,2) ; 
imagesc(squeeze(mean(smthfreqs(:,freqs>40 & freqs<65,:),2))) ; set(gca,'YTick',1:64,'YTickLabel',labs) ;

brwedge = allstims(:,:,1).*allstims(:,:,22) ; 
eptimes = epdata.times(1:10:8700) ; stimtimes = find(eptimes>0 & eptimes <30000) ;  


imagesc(squeeze(smthfreqs(:,15,stimtimes)))
rots = 1:(360/length(stimtimes)):360 ; 
for j=1:10
for i=1:length(rots) ; 
    subplottight(2,2,1) ; 
    imagesc(xg2(:,:,i),[-.5,.5]) ; set(gca,'XTick',[],'YTick',[]) ; 
    subplottight(2,2,3) ; 
    imagesc(xgg2(:,:,i),[-.08,.08]) ; set(gca,'XTick',[],'YTick',[]) ; 
    subplottight(2,2,2) ; imagesc(imrotate(brwedge,-rots(i),'crop')) ; colormap hot ; set(gca,'XTick',[],'YTick',[]) ; 
    getframe 
    
end
end

icount = 1 ;
for i=1:50:length(stimtimes) ; 
    subplottight(3,16,icount) ; 
    imagesc(xg2(:,:,stimtimes(i)),[-.5,.5]) ; set(gca,'XTick',[],'YTick',[]) ; 
    subplottight(3,16,icount+16) ; 
    imagesc(xgg2(:,:,stimtimes(i)),[-.08,.08]) ; set(gca,'XTick',[],'YTick',[]) ; 
    subplottight(3,16,icount+32) ; 
    imagesc(imrotate(brwedge,-rots(i),'crop')) ; colormap hot ; set(gca,'XTick',[],'YTick',[]) ; 
    icount = icount +1 ;  
end









%{
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
allstims(:,:,6:13) = allwedge ; allstims(:,:,14:21) = allbars ; allstims(:,:,22:25) = allquads ; allstims(:,:,26:28) = allrings ; 
%}

