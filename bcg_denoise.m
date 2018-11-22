cd C:\shared\sim_1\first_simultaneous_recordings  ; ls  
clear all ; close all ; 

for dset=1:5 ;
EEG = pop_loadbv('.',['EEG-fMRI_Russell_visual_',num2str(dset),'.vhdr']) ; 
EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
fmri = load_untouch_nii('..\EEG-fMRI 2015-07-30-2\mc_mb4.nii.gz') ; 
boldTR = fmri.hdr.dime.pixdim(5) ; 
eegSR = EEG.srate ; 
tstep = eegSR.*boldTR ; 

taginds = find(strcmp({EEG.urevent.type},'R128')) ; 
lats = cell2mat({EEG.urevent.latency}) ; 
gradstarts = lats(taginds) ; 
clear alldatas ; 
for elec =1:64 ; disp(elec) ; 
        
clear gradepochs meangrads; 
for i=1:length(taginds)-1
   gradepochs(i,:) = EEG.data(elec,gradstarts(i):floor(gradstarts(i)+tstep)-1) ;     
end
meangrads = squeeze(mean(gradepochs,1)) ;  
% subtract out the gradient artifact
new = EEG ; 
for i=1:length(gradstarts)-1 ; 
    new.data(elec,gradstarts(i):floor(gradstarts(i)+tstep)-1) = new.data(elec,gradstarts(i):floor(gradstarts(i)+tstep)-1) - meangrads ; 
end
alldatas(elec,:) = new.data(elec,:) ; 

end
%dend = size(alldatas,2)-5000 ; 
dend = gradstarts(length(gradstarts)) ; 
strt = gradstarts(1)  ;
%strt = 30000 ; 
eeg2 = EEG ; eeg2.data = alldatas ; 
eeg2 = pop_select(eeg2,'point',[strt,dend]) ; 

%%% remove the ballistocardiogram:
eeg2 = pop_resample(eeg2,300) ; 
% get the peaks within a long sliding window(2 seconds) for the correlation
% function
clear maxinds hrepochs hrepochs2
peakdata = eeg2.data(32,:) ; %peakdata = peakdata-smooth(peakdata,2000)' ; 
winlength1 = eeg2.srate*2 ; maxinds = zeros(1,size(eeg2.data,2)) ; 
wprev = floor(eeg2.srate*.33) ; 
wnext = floor(eeg2.srate) ; 
for w=1:winlength1:size(eeg2.data,2)-winlength1
    windowi = peakdata(w:w+winlength1) ; 
    maxi = find(windowi==max(windowi)) ; 
    maxinds(w+maxi) = 1 ; 
end
%plot(eeg2.data(32,:)) ; vline(find(maxinds==1)) ; 
%peakvals = eeg2.data(32,find(maxinds==1)) ; 
peakinds = find(maxinds==1) ; 
for i=2:length(peakinds)-1
    hrepochs(i,:) = peakdata(peakinds(i)-wprev:peakinds(i)+wnext) ; 
end
mhr = squeeze(mean(hrepochs,1)) ; 

%%% get the second peaks (more numerous) and correlate
winlength2 = eeg2.srate ; maxinds = zeros(1,size(eeg2.data,2)) ; 
for w=1:winlength2:size(eeg2.data,2)-winlength2
    windowi = peakdata(w:w+winlength2) ; 
    maxi = find(windowi==max(windowi)) ; 
    maxinds(w+maxi) = 1 ; 
end
peakinds = find(maxinds==1) ; 
for i=1:length(peakinds)-1
    if peakinds(i)-wprev > 0 & peakinds(i) + wnext < length(peakdata)
        hrepochs2(i,:) = peakdata(peakinds(i)-wprev:peakinds(i)+wnext) ; 
    end
end
clear corrs ; 
for i=1:size(hrepochs2,1) ; corrs(i) = corr2(hrepochs2(i,:),mhr) ; end ; 
goods = corrs>0.65 ; 
goodpeakinds = peakinds(goods) ; 
zdiffs = zscore(diff(goodpeakinds)) ; 
tooclose = find(zdiffs<-2) ; goodpeakinds(tooclose) = [] ;
%figure,plot(peakdata) ; vline(goodpeakinds,'r') ; 
geps = hrepochs2(goods,:) ; 
mgep = squeeze(mean(geps(:,wprev-5:wprev+5),2)) ; 
[kc,km] = kmeans(uint8(mat2gray(mgep)*255),4)  ;
for i=1:4 ; kinds{i} = find(km==i) ; meanbcg(i,:) = squeeze(mean(geps(kinds{i},:))) ; end  

eeg3 = eeg2  ; 
clear peakepochs smoothpeaks; 
for i=1:size(eeg2.data,1) ; disp(i) ; 
    for j=2:size(goodpeakinds,2)
        peakepochs(j,:) = eeg2.data(i,goodpeakinds(j)-wprev:(goodpeakinds(j)+wnext)) ;        
    end
    meanpeaks = squeeze(mean(peakepochs,1)) ; 
    for j=2:size(goodpeakinds,2)
        curr = eeg2.data(i,goodpeakinds(j)-wprev:(goodpeakinds(j)+wnext)) ; 
        clear corrs; for c=1:size(peakepochs,1) ; corrs(c) = corr2(curr,peakepochs(c,:)) ; end
        [sv,si] = sort(corrs,'descend') ; closies = si(2:100) ;
        eeg3.data(i,goodpeakinds(j)-wprev:(goodpeakinds(j)+wnext)) = eeg2.data(i,goodpeakinds(j)-wprev:(goodpeakinds(j)+wnext))-mean(peakepochs(closies,:),1) ;    %meanpeaks
        
    end
end

dsets{dset} = eeg3; 
end

% merge the data sets
for i=1:length(dsets) ; 
    if i==1 ; merged = dsets{i} ; 
    elseif i>1 ; merged = pop_mergeset(merged,dsets{i}) ; 
    end
end


%%% POST-PROCESSING 
eeg3 = merged ; 
eeg3.data = eegfiltfft(eeg3.data,eeg3.srate,1,300) ; 
eeg3 = pop_runica(eeg3,'runica') ; 
trigs = {
    'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8','S  9','S 10','S 11','S 12','S 13','S 14','S 15','S 16','S 17','S 18','S 19','S 20'
   %   'S 21','S 22','S 23','S 24','S 25','S 26','S 27','S 28','S 29','S 30','S 31','S 32','S 33','S 34','S 35','S 36','S 37','S 38','S 39','S 40'
} ; 
eps = pop_epoch(eeg3,trigs,[-.8,3.5]) ; 
clear ersp ;
for i=1:size(eps.data,1)
    for j=1:size(eps.data,3)
   [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(eps.icaact(i,:,j)),eps.pnts,[eps.xmin,eps.xmax],eps.srate,0,...
       'plotersp','off','plotitc','off','baseline',NaN,'freqs',[1,120],'nfreqs',60,'winsize',90) ; 
    end
    
end
%figure,for i=1:64 ; subplot(8,8,i) ; topoplot(squeeze(eeg3.icawinv(:,i)),eeg2.chanlocs) ; title(i) ; end
%figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(i,:,:)),[-4,4]) ; title(i) ; end
bersp = zeros(size(ersp)) ; 
for i=1:size(ersp,1)
    for j=1:size(ersp,2)
        for k=1:size(ersp,3)
           mijk = repmat(squeeze(mean(ersp(i,j,k,times<0),4)),[1,size(ersp,4)]) ;
           bersp(i,j,k,:) = squeeze(ersp(i,j,k,:)) - mijk' ; 
        end
    end
end
figure, for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(bersp(i,:,:,:),2)),[-3,3]) ; title(i) ; end

%for i=1:20 ; subplot(4,5,i) ; imagesc(squeeze(ersp(17,i,:,:)),[-40,-10]) ; end



