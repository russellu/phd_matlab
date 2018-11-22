cd c:/shared/denoising_MR ; ls 
clear all ; close all ; 
EEG = pop_loadbv('.','Russell_test_2015-08-05_EEGoutsideMRIroom.vhdr') ; 
EEG = pop_loadbv('.','Russell_test_2015-08-05_EEGnoEPI.vhdr') ; 
EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;

RES = pop_resample(EEG,200) ; 
[t2,f2] = spectopo(RES.data(60,:),0,200,'plot','off') ; 

RES = pop_reref(RES,[]) ; 
%%% remove the ballistocardiogram:
eeg2 = pop_resample(EEG,500) ; 
% get the peaks within a long sliding window(2 seconds) for the correlation
% function
clear maxinds hrepochs hrepochs2
peakdata = eeg2.data(32,:) ;  
winlength1 = eeg2.srate*2 ; maxinds = zeros(1,size(eeg2.data,2)) ; 
wprev = floor(eeg2.srate*.33) ; 
wnext = floor(eeg2.srate) ; 
for w=1:winlength1:size(eeg2.data,2)-winlength1
    windowi = peakdata(w:w+winlength1) ; 
    maxi = find(windowi==max(windowi)) ; 
    maxinds(w+maxi) = 1 ; 
end
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
geps = hrepochs2(goods,:) ; 

eeg3 = eeg2  ; 
clear peakepochs smoothpeaks; 
for i=1:size(eeg2.data,1) ; disp(i) ; 
    if i~=32 ;
    for j=2:size(goodpeakinds,2)
        peakepochs(j,:) = eeg2.data(i,goodpeakinds(j)-wprev:(goodpeakinds(j)+wnext)) ;        
    end
    for j=2:size(goodpeakinds,2)
        curr = eeg2.data(i,goodpeakinds(j)-wprev:(goodpeakinds(j)+wnext)) ; 
        currep = repmat(curr,[size(peakepochs,1),1]) ; 
        ssd = sum((currep-peakepochs).^2,2) ; 
        [sv,si] = sort(ssd,'ascend') ; 
        eeg3.data(i,goodpeakinds(j)-wprev:(goodpeakinds(j)+wnext)) = eeg2.data(i,goodpeakinds(j)-wprev:(goodpeakinds(j)+wnext))-mean(peakepochs(si(2:50),:),1) ;    %meanpeaks      
    end  
    end
    %{
    meanpeaks = squeeze(mean(peakepochs,1)) ; 
    for j=2:size(goodpeakinds,2)
        curr = eeg2.data(i,goodpeakinds(j)-wprev:(goodpeakinds(j)+wnext)) ; 
        clear corrs; for c=1:size(peakepochs,1) ; corrs(c) = corr2(curr,peakepochs(c,:)) ; end
        [sv,si] = sort(corrs,'descend') ; closies = si(2:100) ;
        eeg3.data(i,goodpeakinds(j)-wprev:(goodpeakinds(j)+wnext)) = eeg2.data(i,goodpeakinds(j)-wprev:(goodpeakinds(j)+wnext))-mean(peakepochs(closies,:),1) ;    %meanpeaks      
    end
    %}
end



% post processing
EEG2 = EEG ; 
EEG2 = pop_resample(EEG2,500); 
%EEG2.data = eegfiltfft(EEG2.data,EEG2.srate,1,500) ; 
%EEG2 = pop_runica(EEG2,'runica') ;
clear ersp ; 
trigs = {'S  1','S  4'} ; 
for t=1:length(trigs)
    eps = pop_epoch(EEG2,{trigs{t}},[-.85,2.85]) ;    
    for e=1:size(EEG2.data,1) ;
        for j=1:size(eps.data,3)
            [ersp(t,e,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(eps.data(e,:,j)),eps.pnts,[eps.xmin,eps.xmax],eps.srate,0,...
                'plotersp','off','plotitc','off','baseline',0,'freqs',[1,150],'nfreqs',75,'winsize',130) ; 
        end
    end
end
%figure,for i=1:64 ; subplot(5,13,i) ; topoplot(squeeze(EEG2.icawinv(:,i)),EEG.chanlocs) ; title(i) ; end
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(1,i,:,:)),[-5,5]) ;  title(i) ;end


e1 = squeeze(ersp(1,:,:,:,:)) ; 
for i=1:39 ; subplot(4,10,i) ; imagesc(squeeze(e1(46,i,:,:))) ; end

for i=1:length(goodpeakinds)
   peakepochs(i,:) = eeg2.data(46,goodpeakinds(i):goodpeakinds(i)+eeg2.srate) ; 
    
end
c = corr(peakepochs') ; 
[sv,si] = sort(c(10,:),'descend') ; 


