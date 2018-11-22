clear all ; close all ; 
% problems: need to eliminate the first and last slices (interpolate?)
% => interpolate 2*sliceDurationat the start and end (After aas)

cd c:/shared/badger_eeg/karl ; ls 
vhdrs = dir('*vhdr') ; 
names = {vhdrs.name} ; 
bads = zeros(1,length(vhdrs)) ; 
for i=1:length(names) ; if ~isempty(strfind(names{i},'Pulse')) ; bads(i) = 1 ; end ; end
vhdrs(bads==1) = [] ; 

for v=1:length(vhdrs)
EEG = pop_loadbv('.',vhdrs(v).name) ; 
filt = eegfiltfft(EEG.data(2,:),EEG.srate,50,2500) ; 
mbfactor = 3 ; 
nslices = 33/mbfactor ; 
TR = 0.693 ; 
slicetime = TR/nslices ; 
slicesamples = EEG.srate * slicetime ; 
samplewindow = slicesamples ; 
icount = 1 ; clear epochs inds fullepochs 
allepochs = zeros(ceil(size(EEG.data,2)./samplewindow),size(EEG.data,1),samplewindow+1) ; 
for i=1:samplewindow:size(filt,2)-samplewindow*2
    maxi = find(filt(i:i+samplewindow)==max(filt(i:i+samplewindow))) ; 
    epochs(icount,:) = filt(i+maxi:i+maxi+samplewindow) ; 
    fullepochs(icount,:) = EEG.data(2,i+maxi:i+maxi+samplewindow) ; 
    allepochs(icount,:,:) = EEG.data(:,i+maxi:i+maxi+samplewindow) ; 
    inds(icount) = i+maxi ; 
    icount = icount + 1 ; 
end
mepochs = mean(epochs,1) ; 
corrs = corr(epochs',mepochs') ; 
thresh = 0.95 ; 
goodinds = inds((corrs>thresh)) ; 
tags = zeros(1,size(filt,2)) ; tags(goodinds) = 1 ; 
goodepochs = epochs((corrs>thresh),:) ; 
goodfullepochs = fullepochs(corrs>thresh,:) ; 
allepochs = allepochs(corrs>thresh,:,:) ; 
%{
ssqr = pdist2(goodepochs,goodepochs,'euclidean') ;
navgs=30 ; 
clear topinds
for i=1:size(ssqr,1)
    [sv,si] = sort(ssqr(i,:),'ascend') ; 
    topinds(i,:) = si(2:navgs+1) ; % exclude the top index (itself=0)
end
%}
%{
newchan = EEG.data(2,:) ; nsds = zeros(size(newchan)) ; 
for i=1:length(goodinds)
    artifact = mean(goodfullepochs(topinds(i,:),:),1) ; 
    subbed = EEG.data(2,goodinds(i):goodinds(i)+samplewindow) - artifact ;   
    if i>1
        subdiffs(i) = (subbed(1)-newchan(goodinds(i-1)+samplewindow)) ; 
        nsds(goodinds(i):goodinds(i)+samplewindow) = repmat(subdiffs(i),[numel(artifact),1]) ; 
        subbed = subbed  - (subbed(1)-newchan(goodinds(i-1)+samplewindow)) ; 
    end
    
    prevartifact = artifact ; 
    newchan(goodinds(i):goodinds(i)+samplewindow) = subbed ; 
end
newchan = newchan + smooth(nsds,150)' ; 
interpinds = [goodinds(1)-samplewindow*2:goodinds(1),goodinds(end):goodinds(end)+samplewindow*2] ;
newchan(interpinds) = 0 ; 
%}


inds = 1:size(allepochs,1) ;
newchan = EEG.data ; nsds = zeros(size(newchan)) ;
clear subdiffs 
for i=1:length(goodinds) ; disp(i) ;
    wl = 35 ; 
    artifact = squeeze(mean(allepochs((inds<i+wl & inds>i-wl & inds ~= i),:,:),1)) ; 
    
    epochinds = find(inds<i+wl & inds>i-wl & inds ~= i) ;
    
    filt = repmat(fspecial('gaussian',[1,length(epochinds)],50),[size(allepochs,2),1,size(allepochs,3)]) ; 
    artifact = squeeze(sum(filt.*permute(squeeze(allepochs(epochinds,:,:)),[2,1,3]),2)) ; 
    
    subbed = EEG.data(:,goodinds(i):goodinds(i)+samplewindow) - artifact ;   
    if i>1
        subdiffs(i,:) = (subbed(:,1)-newchan(:,goodinds(i-1)+samplewindow)) ; 
        nsds(:,goodinds(i):goodinds(i)+samplewindow) = repmat(subdiffs(i,:)',[1,size(artifact,2)]) ; 
        subbed = subbed  - repmat(subbed(:,1)-newchan(:,goodinds(i-1)+samplewindow),[1,size(artifact,2)]) ; 
    end
    newchan(:,goodinds(i):goodinds(i)+samplewindow) = subbed ; 
end
newchan = newchan + imfilter(nsds,fspecial('gaussian',[1,150],50)) ; 
interpinds = [goodinds(1)-samplewindow*2:goodinds(1),goodinds(end):goodinds(end)+samplewindow*2] ;
interpinds(interpinds<=0) = [] ; 
newchan(:,interpinds) = 0 ; 
newchan = newchan(:,1:size(EEG.data,2)) ; 

EEG2 = EEG ; 
EEG2.data = newchan ;
EEG2  = pop_resample(EEG2,256) ; raw = EEG2 ; 
pop_saveset(EEG2,['gradient_',strrep(vhdrs(v).name,'.vhdr','.set')]) ; 
[s,f] = spectopo(EEG2.data(:,30000:end-50000),0,EEG2.srate,'plot','off') ; 
figure,imagesc(s) ; 
%EEG2.data = eegfiltfft(raw.data,raw.srate,1,30) ; 
%figure,imagesc(s) ; 

end
%EEG2 = pop_runica(EEG2,'runica') ; 
%raw = applyweights(raw,EEG2) ; 


testchan  = EEG.data(2,:) ; 
filtchan = eegfiltfft(testchan,5000,0,100) ; 




%{
EEG2.data = eegfiltfft(EEG2.data,EEG2.srate,8,16) ; 
EEG2 = pop_runica(EEG2,'runica') ; 


[sica,fica] = spectopo(EEG2.icaact,0,EEG2.srate,'plot','off') ;
EEG2 = pop_chanedit(EEG2,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 


cd outside
outside = pop_loadbv('.','russell_outside.vhdr') ; 
outside = pop_resample(outside,256) ; 
[sout,f] = spectopo(outside.data,0,outside.srate,'plot','off') ; 
subplot(1,2,1) ; imagesc(f,1:64,s,[-40,40]) ; subplot(1,2,2) ; imagesc(f,1:64,sout,[-40,40]) ; 
%}
