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
ecomps = {[14,28,29,45],[4,12,15,34],[12,25,46],[41,42],[30,41],[12,19],[45,39],[31,40,45]} ;

for subby = 1:8 ; 
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

neweeg = pop_loadset('neweeg.set') ; 
merged = ica_applyweights(merged,neweeg) ; 
goods = ecomps{subby} ; bads = zeros(1,64) ; bads(goods) = 1 ; badchans = find(bads==0) ; 
merged = pop_subcomp(merged,badchans) ;
clear ersp ; 
for t=1:length(trigs)
    ep = pop_epoch(merged,trigs{t},[-2,12]) ; 
    for c=1:64 ; 
        for trial=1:size(ep.icaact,3) 
            [ersp(t,c,trial,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.data(c,:,trial)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
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
subplot(3,3,subby) ; 
imagesc(squeeze(mean(mean(anglecomps(ecomps{subby},:,:,:),1),2)))
allcomps(subby,:,:) = squeeze(mean(mean(anglecomps(ecomps{subby},:,:,:),1),2)) ; 

%for i=1:73 ; subplot(8,10,i) ; topoplot(squeeze(mean(mean(anglecomps(:,:,20:30,i),2),3)),merged.chanlocs,'electrodes','off','maplimits',[-55,-20]) ; end
%{
figure,
for i=1:64 ;
    subplot(5,13,i) ; imagesc(imfilter(squeeze(mean(anglecomps(i,:,:,:),2)),fspecial('gaussian',3,3)),[-5,5]) ; title(i);  
end
figure,
for i=1:64 ;
    subplot(5,13,i) ; topoplot(neweeg.icawinv(:,i),neweeg.chanlocs) ; title(i) ; 
end
%}
end

icount = 1 ; clear lowcomps ; 
for i=1:9:73
   lowcomps(:,:,icount) = squeeze(mean(allcomps(:,:,i:i+9.5),3)) ;  
   icount = icount + 1 ; 
end
for s=1:length(subs)
   cd(['c:/shared/anglemasks/',subs{s}]) ;  
   comps = squeeze(lowcomps(s,:,:)) ; 
   save('comps','comps') ;  
end






