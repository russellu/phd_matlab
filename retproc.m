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

subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie','tegan'} ; 
comps= {[13,19,21,22,27,31,34],[7,15,16,17,18,21,22],[7,11,12,17,19,30],[7,23],[8,16,21],[20,27],[11,25]} ; 
for subby =8%:7 ; 
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
sounds=dir('highfreq*allstim*set') ;
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

figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(anglecomps(i,:,:,:),2)),[-5,5]) ; title(i);  end
%{
figure('Position',[100,100,1500,400]) ; 
for i=1:length(comps{subby})
   subplot(3,8,i) ; 
   imagesc(uniques,freqs,squeeze(mean(anglecomps(comps{subby}(i),:,:,:),2)),[-3,3]) ;  axis xy ; xlabel('angle(deg)') ;ylabel('frequency(hz)') ; 
   subplot(3,8,i+8) ; 
   topoplot(ep.icawinv(:,comps{subby}(i)),ep.chanlocs) ; 
   subplot(3,8,i+16) ; 
   plot(squeeze(mean(mean(anglecomps(comps{subby}(i),:,freqs>8 & freqs<20,:),2),3))) ; hold on ; plot(squeeze(mean(mean(anglecomps(comps{subby}(i),:,freqs>40 & freqs<60,:),2),3)),'r') ;
   xlim([0,size(anglecomps,4)]) ; hline(0,'k') ;  set(gca,'XTick',1:25:length(uniques),'XTickLabel',uniques(1:25:end)) ; xlabel('angle(deg)') ; ylabel('power(db)') ;
end
%}
end

uniques = imresize(uniques,[size(anglecomps,4),1]) ; 
figure, for i=1:length(comps{subby}) ; subplot(2,4,i) ; imagesc(abs(squeeze(mean(anglecomps(comps{subby}(i),:,:,:),2))),[-4,4]); end

cd(['c:/shared/badger_mri/',subs{subby},'/nii']) ; ls  
fmri = load_untouch_nii('allangles.nii.gz') ; fimg = fmri.img ; 
manglecomps = squeeze(mean(anglecomps,2)) ; 
clear resangles
for i=1:64 ; resangles(i,:,:) = imresize(squeeze(manglecomps(i,:,:)),[60,size(fimg,4)]) ; end
newcomps = {[2,6],[4,1],[1,6],[1,2],[1,2],[2,1],[1,2]} ; 







[p,s] = polyfit(1:90,squeeze(mean(abs(resangles(20,3:30,:)),2))',5) ;
x=1:90 ; 
y1 = polyval(p,x) ; plot(y1) ; 
figure,subplot(1,2,1) ; plot(squeeze(mean(abs(resangles(20,3:30,:)),2))') ; subplot(1,2,2); plot(y1) ; 
% do correlation stuff:


cd(['c:/shared/badger_mri/',subs{subby},'/nii']) ; ls  
% get the fmri correlation mask 
anat = load_untouch_nii('f_topup_mc_retino_allstims_01.nii.gz') ; anatim = anat.img ;
delayavg = load_untouch_nii('delay_mean.nii.gz') ; 
stdavg = squeeze(delayavg.img(:,:,:,1,2)) ; 
[kc,km] = kmeans(uint8(mat2gray(reshape(stdavg,[size(stdavg,1),size(stdavg,2)*size(stdavg,3)]))*255),2) ; 
masky = reshape(km==2,size(stdavg)) ; 
anat.img = masky ; save_untouch_nii(anat,'masky.nii.gz') ; 



figure,
clear corrs 
for c=1:length(newcomps{subby})
    
   [p,s] = polyfit(1:90,squeeze(mean(abs(resangles(comps{subby}(newcomps{subby}(c)),3:30,:)),2))',5) ;
    x=1:90 ; 
    y1 = polyval(p,x) ; subplot(2,4,c) ; plot(y1) ; 
   corrs(c,:,:,:) = voxcorr(fimg,y1) ;  
    
end
clear corrrgb 
corrrgb(:,:,:,1) = squeeze(corrs(1,:,:,:)).*masky ; 
corrrgb(:,:,:,3) = squeeze(corrs(2,:,:,:)).*masky ; 
corrrgb(isnan(corrrgb)) = 0 ;corrrgb(isinf(corrrgb)) = 0 ; 
corrrgb(corrrgb<0) = 0  ; corrrgb = mat2gray(corrrgb) ; 

%{
for c=1:length(comps{subby}) ; disp(c) ; 
    for f=1:60
        corrs(c,f,:,:,:) = abs(voxcorr(fimg,squeeze(resangles(comps{subby}(c),f,:)))) ; 
    end
end
corrlow = (squeeze(mean(corrs(:,freqs>8 & freqs<18,:,:,:),2))) ; corrlow = permute(corrlow,[2,3,4,1]) ; corrlow(isnan(corrlow)) = 0 ; 
corrhigh = (squeeze(mean(corrs(:,freqs>40 & freqs<70,:,:,:),2))) ; corrhigh = permute(corrhigh,[2,3,4,1]) ; corrhigh(isnan(corrhigh)) = 0 ; 
clear corrrgb ;
corrrgb(:,:,:,1) = squeeze(corrlow(:,:,:,newcomps{subby}(1))).*double(squeeze(corrlow(:,:,:,newcomps{subby}(1))>0)) + squeeze(corrhigh(:,:,:,newcomps{subby}(1))).*double(squeeze(corrhigh(:,:,:,newcomps{subby}(1))>0)) ; 
corrrgb(:,:,:,3) =  squeeze(corrlow(:,:,:,newcomps{subby}(2))).*double(squeeze(corrlow(:,:,:,newcomps{subby}(2))>0)) + squeeze(corrhigh(:,:,:,newcomps{subby}(2))).*double(squeeze(corrhigh(:,:,:,newcomps{subby}(2))>0)) ; 
corrboth = (corrlow + corrhigh) / 2 ; 
%}
figure,for i=1:length(comps{subby}) ; subplot(2,4,i) ; topoplot(squeeze(ep.icawinv(:,comps{subby}(i))),ep.chanlocs) ; end 
anat = load_untouch_nii('f_topup_mc_retino_allstims_01.nii.gz') ; anatim = anat.img ;
%for zs=10:20 ; figure,
%for i=1:length(comps{subby}) ; subplot(2,4,i) ; plotoverlayIntensity2D(squeeze(mean(anatim(:,:,zs),3)),(abs(squeeze(mean(corrboth(:,:,zs,i),3)))),squeeze(mean(corrboth(:,:,zs,i),3)),270) ; end
%end


%anat.img = squeeze(corrrgb(:,:,:,1)) ; save_untouch_nii(anat,'comp1.nii.gz') ; 
%anat.img = squeeze(corrrgb(:,:,:,3)) ; save_untouch_nii(anat,'comp2.nii.gz') ; 

% make the retinotopic tuning
for i=1:64; for j=1:60 ; manglecomps(i,j,:) = smooth(squeeze(manglecomps(i,j,:))) ; end ; end

% ground truth:






[xg,yg] = meshgrid(-200:199,-200:199) ; 
[theta,rho] = cart2pol(xg,yg) ; mask = rho<200 ; 
degangles = (flipud(mod((theta*180)/pi+360,360))) ; 
plotfreqs = {freqs>10 & freqs<25, freqs>40 & freqs<60} ; 
titles1 = {'comp1 alpha/beta tuning', 'comp2 alpha/beta tuning'} ;
titles2 = {'comp1 gamma tuning', 'comp2 gamma tuning'} ;
for pf = 1:length(plotfreqs) ; 
    gfs = plotfreqs{pf} ; 
    zangles = zeros(size(degangles)) ; 
    for c=1:length(newcomps{1})
        for i=1:length(uniques)
            if i==1 ;
                zangles(degangles<uniques(i+1) & degangles>=0) = squeeze(abs(mean(manglecomps(comps{subby}(newcomps{subby}(c)),gfs,i),2))) ; 
            else
                zangles(degangles<=uniques(i) & degangles>uniques(i-1)) = squeeze(abs(mean(manglecomps(comps{subby}(newcomps{subby}(c)),gfs,i),2))) ; 
            end
        end
        zangles(degangles>=uniques(length(uniques))) = (squeeze(abs(mean(manglecomps(comps{subby}(newcomps{subby}(c)),gfs,i),2))) + squeeze(abs(mean(manglecomps(comps{subby}(newcomps{subby}(c)),gfs,1),2))))/2 ;
        if pf==1 ; subplot(2,2,c) ; imagesc(zangles) ; colormap hot ; set(gca,'XTickLabel',[],'YTickLabel',[]) ; title(titles1{c}) ; 
        else subplot(2,2,c+2) ; imagesc(zangles) ; colormap hot ; set(gca,'XTickLabel',[],'YTickLabel',[]) ; title(titles2{c}) ; 
        end   
       allzangles(pf,c,:,:) = zangles ; 
    end
end

clear rgbz
rgbz(:,:,1) = squeeze(allzangles(1,1,:,:) + allzangles(2,1,:,:)) ; 
rgbz(:,:,3) = squeeze(allzangles(1,2,:,:) + allzangles(2,2,:,:)) ; 
figure,imagesc(uint8(mat2gray(rgbz)*255)) ; 






for i=1:64 ; [~,g(i,:,:)] = topoplot(EEG.icawinv(:,i),EEG.chanlocs) ; end
fhandle = figure('Position',[100,100,1000,800]) ; 
subplot(2,2,1) ; 
imagesc(uniques,freqs,abs(squeeze(manglecomps(comps{subby}(newcomps{subby}(1)),:,:)))) ;  colorbar ; colormap jet ; axis xy ; xlabel('polar angle(rad)') ; ylabel('frequency(hz)') ;  title('comp1 retinotopic tuning') ; 
subplot(2,2,2) ; 
imagesc(uniques,freqs,abs(squeeze(manglecomps(comps{subby}(newcomps{subby}(2)),:,:)))) ; colorbar ; colormap jet ; axis xy ;   title('comp2 retinotopic tuning') ; 
clear rgb
rgb(:,:,1) = imfilter((abs(squeeze(manglecomps(comps{subby}(newcomps{subby}(1)),:,:)))),fspecial('gaussian',3,3)) ;
rgb(:,:,3) = imfilter((abs(squeeze(manglecomps(comps{subby}(newcomps{subby}(2)),:,:)))),fspecial('gaussian',3,3)) ;
rgb = uint8(mat2gray(rgb)*255) ; 
subplot(2,2,3) ; 
%imagesc(uniques,freqs,rgb) ; axis xy ; xlabel('polar angle(deg)') ; ylabel('frequency(hz)') ;  title('comp1=red, comp2=blue, retinotopic tuning') ; 
imagesc(uint8(mat2gray(rgbz)*255)) ; title('red=comp1,blue=comp2, modulation(polar angle)') ; set(gca,'XTick',[],'YTick',[]) ; 
subplot(2,2,4) ; 
clear rgbtopo
rgbtopo(:,:,1) = squeeze(flipud(squeeze(g(comps{subby}(newcomps{subby}(1)),:,:)))) ; 
rgbtopo(:,:,3) = squeeze(flipud(squeeze(g(comps{subby}(newcomps{subby}(2)),:,:)))) ; 
rgbtopo(isnan(rgbtopo)) = 1 ; 
imagesc(uint8(mat2gray(rgbtopo)*255))  ; title('scalp topography, comp1=red, comp2=blue') ;  set(gca,'XTick',[],'YTick',[]) ; 

figure,
for i=1:15; subplot(3,5,i)  ;
plotbrainRGB(squeeze(anatim(:,:,i+3)),mat2gray(squeeze(mean(corrrgb(:,:,i+3,:),4))),squeeze(corrrgb(:,:,i+3,:)),270) ;
end


%{
cd c:/shared/badger_mri/dina/nii/  ; t1 = load_untouch_nii('t1.nii') ; 
t1im = t1.img ;icount=1 ;
for i=1:10:150; subplot(3,5,icount)  ;
plotoverlayIntensity2D(squeeze(t1im(:,:,i+3)),mat2gray(squeeze(t1im(:,:,i+3))),squeeze(t1im(:,:,i+3)),270) ;
icount=icount+1;
end
%}