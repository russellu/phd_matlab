clear all ; close all 
name = 'alex' ; 

cd(['c:/shared/badger_eeg/',name]) ; ls 

% get a raw epoch for the times
setnames = dir('preproc*gamma*set') ; eegsets = {setnames(1).name,setnames(2).name} ; EEG = pop_loadset(eegsets{1}) ; ep = pop_epoch(EEG,{'S  1'},[-6,12]) ; 
baseline = find(ep.times<0 & ep.times>-3000) ; task = find(ep.times>0 & ep.times<5000) ; 

% try convolving the FMRi with and without the baseline EEG power..

% load and visualize the EEG single trial data
eeg = load('allepochs.mat') ; eeg = eeg.allepochs ; 
meanfreq = squeeze(mean(mean(eeg(1,:,:,:,:),3),5)) ; 
meanelecs = squeeze(mean(eeg,3)) ; % power
meanelecs = log(meanelecs.^2) ; % squared power
% mean squared power
meantask = squeeze(mean(mean(eeg(:,:,:,task,:).^2,3),4)) ; 
meanbase = squeeze(mean(mean(eeg(:,:,:,baseline,:).^2,3),4)) ; 
%icount=1 ; for i=1:5:120 ; subplot(3,8,icount) ; bar(squeeze(mean(meanbase(:,i:i+5,:),2))) ; icount = icount +1 ; title(i) ; end
hrf = spm_hrf(0.693) ;

resmean = zeros(size(meanelecs,1),size(meanelecs,2),ceil(18/0.693),size(meanelecs,4)) ; 
conved = zeros(size(resmean)) ; 
% resize the EEG to the FMRI sampling rate:
for i=1:size(meanelecs,1) ; disp(i) ; 
    for j=1:size(meanelecs,2)
        for k=1:size(meanelecs,4)
            x = imresize(squeeze(meanelecs(i,j,:,k)),[ceil(18/0.693),1]) ; 
            resmean(i,j,:,k) = x ; 
            convx = conv(x,hrf,'full') ; 
            convx = convx(1:length(x)) ; % resize to original
            conved(i,j,:,k) = convx ;            
        end
    end
end


cd(['C:\shared\badger_mri\',name,'\nii']) ; 
allts = load('allts.mat') ; allts = allts.allts ; 
peakinds = load('peakinds.mat') ; peakinds = peakinds.peakinds ; taskf = load('task.mat') ; taskf = taskf.task ; 

res = permute(allts,[1,6,3,4,5,2]) ; 
res96 = zeros(18,size(res,3),size(res,4),size(res,5),96) ;
res96(:,:,:,:,1:32) = squeeze(res(1,:,:,:,:,:)) ; res96(:,:,:,:,33:64) = squeeze(res(1,:,:,:,:,:)) ; res96(:,:,:,:,65:96) = squeeze(res(1,:,:,:,:,:)) ; 

shortconved = conved(:,:,size(conved,3)-17:end,:) ; 
resconv = zeros(120,18,96) ; 
resconv(:,:,1:32) = squeeze(shortconved(1,:,:,:)) ; resconv(:,:,33:64) = squeeze(shortconved(2,:,:,:)) ; resconv(:,:,65:96) = squeeze(shortconved(3,:,:,:)) ; 

for xx=1:15 ; 
clear corrs ; 
icount = 1 ; 
for i=1:5:120-5
    corrs(:,:,:,icount) = voxcorr(squeeze(mean(res96(10:14,:,:,:,:),1))-squeeze(mean(res96(1:3,:,:,:,:),1)),squeeze(mean(mean(resconv(i:i+5,xx:xx+1,:),1),2))) ; 
    icount = icount + 1 ; 
end
f = load_untouch_nii('f_retino_allstims_01.nii') ; 
zinds = 25:35; figure,
for i=1:size(corrs,4) ; subplot(5,5,i),
    plotoverlayIntensity2D(squeeze(mean(f.img(zinds,:,:),1)),mat2gray(abs(squeeze(mean(corrs(zinds,:,:,i),1)))),squeeze(mean(corrs(zinds,:,:,i),1)),90) ; title(i*5) ; 
end
end


%{
% get the bad indices in the EEG data:
clear goods
for i=1:3 ;
    for j=1:120
        zij = squeeze(zscore(mean(shortconved(i,j,10:14,:),3))) ; 
        goods{i,j} = find(abs(zij)<2) ; 
    end
end

% for the convolved power
clear corrbrains
hcount=1 ;
for hz=1:5:120-5 ; 
    for i=1:3 ;
        corrbrains(i,hcount,:,:,:) = voxcorr((squeeze(mean(res(i,10:14,:,:,:,:),2))-squeeze(mean(res(i,1:3,:,:,:,:),2))),...
            squeeze(mean(mean(shortconved(i,hz:hz+5,peakinds,:),3),2))-squeeze(mean(mean(shortconved(i,hz:hz+5,1:3,:),3),2))) ; 
    end
    hcount = hcount +1 ; 
end
%for i=1:5:120-5 ; figure ; imagesc(squeeze(mean(mean(mean(corrbrains(3,i:i+5,:,:,8:15),1),5),2)),[-.5,.5]) ; title(i) ; end

f = load_untouch_nii('f_retino_allstims_01.nii') ; 

zinds = 25:35; figure,
for i=1:size(corrbrains,2) ; subplot(5,5,i),
    plotoverlayIntensity2D(squeeze(mean(f.img(zinds,:,:),1)),mat2gray(abs(squeeze(mean(corrbrains(3,i,zinds,:,:),3)))),squeeze(mean(corrbrains(3,i,zinds,:,:),3)),90) ; title(i*5) ; 
end

%}
%{
mkdir corrs ; cd corrs ; 
for i=1:5:120-5 ;  
    nii = (squeeze(mean(mean(mean(corrbrains(1:3,i:i+5,:,:,:),1)),2))) ;
    f.img = single(nii) ;   
    save_untouch_nii(f,[num2str(i),'.nii.gz']) ; 
end
%}
% something else to try: convolve the stimulus time course in each
% frequency band for each trial and compare it to BOLD. 













%{
% for the raw power values:
clear corrbrains
for hz=1:120 ; 
    for i=1:3 ;
        corrbrains(i,hz,:,:,:) = voxcorr((squeeze(mean(res(i,10:14,:,:,:,:),2))-squeeze(mean(res(i,1:3,:,:,:,:),2))),...
                                        ((squeeze(mean(meantask(i,hz,:),2)))-(squeeze(mean(meanbase(i,hz,:),2)))) ./ (squeeze(mean(meanbase(i,hz,:),2))) ) ;   % 
    end
end
for i=1:5:120-5 ; figure ; imagesc(squeeze(mean(mean(mean(corrbrains(3,i:i+5,:,:,5:10),1),5),2)),[-.5,.5]) ; title(i) ; end
%}
%{
zs = zscore(sum(meanfreq,2)) ; goodfreq = find(abs(zs)<8) ; 
for i=3:size(meanfreq,1)
   plot(meanfreq(i,:),'Color',[i/120,0,0]) ;  
   hold on 
end
for trial=1:16
    subplot(4,4,trial) ; 
    for i=3:size(meanfreq,1)
        plot(squeeze(mean(allepochs(1,i,:,:,trial),3)),'Color',[i/120,0,0]) ; 
       hold on 
    end
end
for i=1:3 ; figure ; for j=1:32 ; subplot(4,8,j) ; imagesc(log(squeeze(mean(eeg(i,:,:,:,j),3)).^2)) ; end ; end
figure,
subplot(1,2,1) ; imagesc(squeeze(std(log(meantask),0,3))) ; title('standard deviation task') ; 
subplot(1,2,2) ; imagesc(squeeze(std(log(meanbase),0,3))) ; title('standard deviation base') ; 
stdlogbase = squeeze(std(log(meanbase),0,3)) ; 
stdlogtask = squeeze(std(log(meantask),0,3)) ; 
clear corrs ;
for i=1:3
    for j=1:120 
        corrs(i,j,:,:) = corr([squeeze(meanbase(i,j,:)),squeeze(meantask(i,j,:))]) ; 
    end
end
%}

