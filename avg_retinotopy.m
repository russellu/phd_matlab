%%%% process the EEG retinotopic mapping
% quadrant starts in the bottom left, and is rotated by startangles
clear all ; close all ;
% so stimuli are the following:
trigs{1} = {'S 11','S 12','S 13','S 14'} ; 
trigs{2} = {'S 21','S 22','S 23','S 24'} ; 
trigs{3} = {'S 31','S 32','S 33','S 34'} ; 
trigs{4} = {'S 41','S 42','S 43','S 44'} ; 
trigs{5} = {'S 51','S 52','S 53','S 54'} ; 
trigs{6} = {'S 61','S 62','S 63','S 64'} ; 

subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie'} ; 
for subby = 1:6%7 ; 
startangles = [0,60,120,180,240,300] ;

if subby==1 || subby==4
    realangles = mod(225-startangles,360) ; 
    r = ones(1,6) ; 
    radangles = (realangles*pi)/180 ; 
    [x,y] = pol2cart(radangles,r) ; 
    [theta,rho] = cart2pol(-x,y) ; 
    deg = theta*180/pi ; 
    realangles = mod(deg,360) ;    
else
    realangles = mod(225-startangles,360) ; 
end
cd(['c:/shared/badger_eeg/',subs{subby}]) ; ls  ;
sounds=dir('1Hz*allstim*set') ;
for i=1:max(size(sounds)) ;  
   EEG = pop_loadset(sounds(i).name) ; 
   if i==1
      merged = EEG ; 
   else merged = pop_mergeset(EEG,merged,1) ; 
   end
end

clear ersp ; 
for t=1:length(trigs)
    ep = pop_epoch(merged,trigs{t},[-2,12]) ; 
    for c=1:64 ; 
        for trial=1:size(ep.icaact,3) 
            [ersp(t,c,trial,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(c,:,trial)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                        'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'timesout',100) ; 
        end
    end
end

bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,size(ersp,5)]) ; 
mbersp = squeeze(mean(bersp,3)) ; 
%for i=1:6 ; figure ; for j=1:64 ; subplot(5,13,j) ; imagesc(squeeze(mean(bersp(i,j,:,:,:),3)),[-8,8]) ; end ; end 
stiminds = find(times>0 & times<10) ; % indices during which the stimulus was present
nangles = length(stiminds) ; 
% start angle (0) is taken as 180 + 45 = 225
angleincr = 360/nangles ; 
for i=1:length(realangles)
    if subby==1 || subby==4
        fullangles(i,:) = mod(realangles(i):angleincr:realangles(i)+360,360) ;
    else
        fullangles(i,:) = mod(realangles(i)+360:-angleincr:realangles(i),360) ; % negative increment for clockwise rotation
    end
end
uniques = unique(fullangles) ; 
anglecomps = zeros(size(bersp,2),size(bersp,3),size(bersp,4),size(fullangles,2)) ; % components, freqs, angles
for i=1:size(bersp,1) % for all starting angles
    [~,si] = sort(fullangles(i,:),'descend') ; 
    si = si + sum(times<0) ; % adjust indices for baseline
    % si is now the angle indices for power at that angle
    for j=1:size(bersp,2) % for all components
        berspij = squeeze(bersp(i,j,:,:,si)) ; 
        anglecomps(j,:,:,:) = squeeze(anglecomps(j,:,:,:)) + berspij./6 ; 
    end
end
allanglecomps(subby,:,:,:,:) = anglecomps ; 
eegs{subby} = EEG ; 
%figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(anglecomps(i,:,:,:),2)),[-3,3]) ; title(i);  end
%tp(EEG) ; 
end
uniques = imresize(uniques,[size(anglecomps,4),1]) ; 

comps = {[31,19],[16,17],[30,7],[7,23],[16,8],[20,27]} ; 

for i=1:length(comps)
    subangles(i,1,:,:,:) = squeeze(allanglecomps(i,comps{i}(1),:,:,:)) ; 
    subangles(i,2,:,:,:) = squeeze(allanglecomps(i,comps{i}(2),:,:,:)) ; 
end
msubs=  squeeze(mean(subangles,1)) ; 

subplot(2,2,1) ;  imagesc(squeeze(mean(msubs(1,:,:,:),2)),[-3,3])
subplot(2,2,2) ;  imagesc(squeeze(mean(msubs(2,:,:,:),2)),[-3,3])
subplot(2,2,3) ; 
plot(squeeze(mean(mean(msubs(1,:,freqs>8 & freqs<15,:),2),3))) ; hold on ; hline(0,'k') ; 
plot(squeeze(mean(mean(msubs(1,:,freqs>35 & freqs<65,:),2),3)),'r') ; xlim([0,74]) ; ylim([-3.5,1.5]) ;
subplot(2,2,4) ; 
plot(squeeze(mean(mean(msubs(2,:,freqs>8 & freqs<15,:),2),3))) ; hold on ; hline(0,'k') ; 
plot(squeeze(mean(mean(msubs(2,:,freqs>35 & freqs<65,:),2),3)),'r') ; xlim([0,74]) ; ylim([-3.5,1.5]) ;

mmsubs = squeeze(mean(mean(msubs))) ; 
figure,subplot(1,2,1) ; imagesc(mmsubs,[-3,3]) ;
subplot(1,2,2) ; plot(mean(squeeze(mmsubs(freqs>8 & freqs<18,:)))) ; hold on ; 
plot(mean(squeeze(mmsubs(freqs>35 & freqs<60,:))),'r') ; hline(0,'k') ; 

mtsubs = squeeze(mean(msubs,2)) ; 


for c=1:2 
[xg,yg] = meshgrid(-200:199,-200:199) ; 
[theta,rho] = cart2pol(xg,yg) ; mask = rho<200 ; 
degangles = (flipud(mod((theta*180)/pi+360,360))) ; 
plotfreqs = {freqs>10 & freqs<25, freqs>40 & freqs<60} ; 
for pf = 1:length(plotfreqs) ; 
    gfs = plotfreqs{pf} ; 
    zangles = zeros(size(degangles)) ; 
    for i=1:length(uniques)
        if i==1 ;
            zangles(degangles<uniques(i+1) & degangles>=0) = squeeze(mean(mtsubs(c,plotfreqs{pf},i))) ; 
        else
            zangles(degangles<=uniques(i) & degangles>uniques(i-1)) = squeeze(mean(mtsubs(c,plotfreqs{pf},i))) ; 
        end
    end
    zangles(degangles>=uniques(length(uniques))) = (squeeze(mean(mtsubs(c,plotfreqs{pf},i))) + mean(squeeze(mtsubs(c,plotfreqs{pf},1))))/2 ; 
   allzangles(pf,:,:) = zangles ; 
end
figure,
subplot(1,2,1) ; imagesc(abs(squeeze(allzangles(1,:,:))),[0,3]) ; colormap hot
subplot(1,2,2) ; imagesc(abs(squeeze(allzangles(2,:,:))),[0,3]) ; colormap hot

zvals(c,1,:,:) = abs(squeeze(allzangles(1,:,:))) ; 
zvals(c,2,:,:) = abs(squeeze(allzangles(2,:,:))) ; 

end

rgb1(:,:,1) = squeeze(zvals(1,1,:,:)) ;
rgb1(:,:,3) = squeeze(zvals(2,1,:,:)) ; 

rgb2(:,:,1) = squeeze(zvals(1,2,:,:)) ;
rgb2(:,:,3) = squeeze(zvals(2,2,:,:)) ; 

subplot(1,2,1) ; imagesc(uint8(mat2gray(rgb1)*255)) ; title('alpha/beta') ; 
subplot(1,2,2) ; imagesc(uint8(mat2gray(rgb2)*255)) ; title('gamma') ;

clear icawinvs
for i=1:length(eegs)
   icawinvs(1,i,:) = squeeze(eegs{i}.icawinv(:,comps{i}(1))) ;  
   icawinvs(2,i,:) = squeeze(eegs{i}.icawinv(:,comps{i}(2))) ;      
end

subplot(1,2,1) ; [~,g(1,:,:)] = topoplot(squeeze(mean(icawinvs(1,:,:),2)),EEG.chanlocs) ; 
subplot(1,2,2) ; [~,g(2,:,:)] = topoplot(squeeze(mean(icawinvs(2,:,:),2)),EEG.chanlocs) ; 

rgbtopo(:,:,1) = flipud(squeeze(g(1,:,:))) ; 
rgbtopo(:,:,3) = flipud(squeeze(g(2,:,:))) ; rgbtopo(isnan(rgbtopo)) = 0 ; 
imagesc(uint8(mat2gray(rgbtopo)*255))




