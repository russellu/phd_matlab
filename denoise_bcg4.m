function EEG2 = denoise_bcg4(EEG) 
% function EEG2 = denoise_bcg4(EEG) ;

EEGfilt = EEG ; EEGfilt.data = eegfiltfft(EEG.data,EEG.srate,5,15) ; 
corrmat = corr(EEGfilt.data') ; 
stelec = 46 ; 
[~,si] = sort(corrmat(stelec,:),'descend') ;
meanpost = (mean(EEGfilt.data(si(1:10),:),1)) ; wsize2 = round(EEG.srate*.25) ;

clear diffepochs
for elec =1:64 ;
wsize = round(EEG.srate*1.5) ; icount = 1 ; clear maxinds mepochs; 
for i=wsize:wsize:size(EEGfilt.data,2)-wsize*2
    maxinds(icount) = find(meanpost(i:i+wsize)==max(meanpost(i:i+wsize)))+i ; 
    mepochs(icount,:) = meanpost(maxinds(icount)-round(wsize2):maxinds(icount)+round(wsize2)) ; 
    icount = icount + 1 ; 
end
meanepoch = mean(mepochs,1) ; 

icount = 1 ; clear maxinds mepochs corrmepochs ; 
for i=wsize2:wsize2:size(EEGfilt.data,2)-wsize2*8
    maxinds(icount) = find(meanpost(i:i+wsize)==max(meanpost(i:i+wsize)))+i ; 
    mepochs(icount,:) = EEG.data(elec,maxinds(icount)-round(wsize/4):maxinds(icount)+round(wsize/2)) ; 
    corrmepochs(icount,:) = meanpost(maxinds(icount)-round(wsize2):maxinds(icount)+round(wsize2)) ; 
    icount = icount + 1 ; 
end
mepochs(diff(maxinds)==0,:) = [] ; corrmepochs(diff(maxinds)==0,:) = [] ; maxinds(diff(maxinds)==0) = []  ;
meancorrs = corr(corrmepochs',meanepoch') ; 
badzs = find(zscore(meancorrs)<-1) ; 
maxinds(badzs) = [] ; 
mepochs(badzs,:) = [] ; 
plot(EEG.data(46,:)) ; vline(maxinds) ;

avamt = 20 ;
for ep=1:size(mepochs,1) ; 
    cmat = corr(mepochs') ; 
    [~,si] = sort(cmat(ep,:),'descend') ; 
    diffepochs(elec,ep,:) = mepochs(si(1),:) - mean(mepochs(si(2:avamt),:)) ; 
end

disp(elec) ; 
end

% replace the indices
EEG2 = EEG ; 
for i=1:length(maxinds)
   EEG2.data(:,maxinds(i)-round(wsize/4):maxinds(i)+round(wsize/2)) = squeeze(diffepochs(:,i,:)) ; 
end

end