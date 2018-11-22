%cd('c:/Vision/Raw Files') ; ls 
cd('c:/shared/denoising_MR') ; 
clear all ; close all ; 
EEGnoise = pop_loadbv('.','Russell_test_2015-08-05_EEGnoEPI.vhdr') ; 
EEGnoise = pop_chanedit(EEGnoise,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
eeg2 = EEGnoise ; eeg2 = pop_resample(eeg2,500) ; 
srate = eeg2.srate ; 
%eeg2.data = eegfiltfft(eeg2.data,eeg2.srate,1,250) ; 
% get the peaks within a long sliding window(2 seconds) for the correlation
% function
clear maxinds hrepochs hrepochs2
peakdata = eeg2.data(32,:) ;  
winlength1 = eeg2.srate*2 ; maxinds = zeros(1,size(eeg2.data,2)) ; 
wprev = 100*(500/srate) ; 
wnext = 500*(500/srate) ; 
for w=1:winlength1:size(eeg2.data,2)-winlength1
    windowi = peakdata(w:w+winlength1) ; 
    maxi = find(windowi==max(windowi)) ; 
    maxinds(w+maxi) = 1 ; 
end
hpeakinds = find(maxinds==1) ; 
for i=2:length(hpeakinds)-1
    hrepochs(i,:) = peakdata(hpeakinds(i)-wprev:hpeakinds(i)+wnext) ; 
end
mhr = squeeze(mean(hrepochs,1)) ; 

%%% get the second peaks (more numerous) and correlate
winlength2 = eeg2.srate ; maxinds = zeros(1,size(eeg2.data,2)) ; 
for w=1:round(winlength2):size(eeg2.data,2)-winlength2
    windowi = peakdata(w:w+winlength2) ; 
    maxi = find(windowi==max(windowi)) ; 
    maxinds(w+maxi) = 1 ; 
end
hpeakinds = find(maxinds==1) ; 
for i=1:length(hpeakinds)-1
    if hpeakinds(i)-wprev > 0 & hpeakinds(i) + wnext < length(peakdata)
        hrepochs2(i,:) = peakdata(hpeakinds(i)-wprev:hpeakinds(i)+wnext) ; 
    end
end
clear corrs ; 
for i=1:size(hrepochs2,1) ; corrs(i) = corr2(hrepochs2(i,:),mhr) ; end ; 
goods = corrs>0.65 ; goodepochs = hrepochs2(goods,:) ; 
goodpeakinds = hpeakinds(goods) ; 
% find the double peak indices and remove the lower one
zdiffs = zscore(diff(goodpeakinds)) ; 
tooclose = find(zdiffs<-2) ; 
for i=1:length(tooclose)
    if peakdata(goodpeakinds(tooclose(i))) < peakdata(goodpeakinds(tooclose(i))-1)       
        tooclose(i) = (tooclose(i)+1) ;  
    end
end
goodepochs(tooclose,:) = [] ; goodpeakinds(tooclose) = [] ; 
geps = goodepochs ; 
%plot(eeg2.data(32,:)) ; vline(goodpeakinds)
imagesc(geps);

for i=1:length(goodpeakinds)
    scalpepochs(i,:,:) = eeg2.data(:,goodpeakinds(i)-wprev:goodpeakinds(i)+wnext) ; 
     
end

scalp = squeeze(mean(scalpepochs,1)) ; 
es = elocs ; 
est = zeros(size(es,1),size(es,2),size(scalp,2)) ; 
for i=1:size(es,1)
    for j=1:size(es,2)
        if es(i,j) == 0 ; est(i,j,:) = zeros(1,size(scalp,2)) ;
        else est(i,j,:) = scalp(es(i,j),:) ; 
        end    
    end
end
disp3d(est) ; 

%x = zeros(1,size(scalpepochs,2),size(scalpepochs,3)) ; 
%x(1,:,:) = scalp ; 
%save_nii(make_nii(scalpepochs),'scalpepochs.nii.gz') ; 
%save_nii(make_nii(repmat(x,[size(scalpepochs,1),1,1])),'meanscalpepochs.nii.gz') ; 



for i=1:size(scalpepochs,3)
   scorrs(:,:,i) = corr(squeeze(scalpepochs(:,:,i))) ;  
end
scorrs(scorrs==1) = 0 ; 
% reconstruct the heart beat based on a linear combination of all other
% channels

for trial=1:size(scalpepochs,1)
    for tp=1:size(scalpepochs,3)
        etrials(trial,tp) = squeeze(scalpepochs(trial,:,tp))*squeeze(scorrs(1,:,tp))' ; 
    
    
    end
end


 
% subtract out the n-highest correlating channels at each time point.
for i=1:size(scorrs,3)
    for j=1:64
        [sv(i,j,:),si(i,j,:)] = sort(squeeze(scorrs(:,j,i)),'descend') ; 
    end
end
for i=1:size(scalpepochs,1) 
    for j=1:size(scalpepochs,2)
        subepochs(i,j,:) = squeeze(scalpepochs(i,j,:))-squeeze(mean(sv(:,j,1:10),3)) ; 
        
    end
end


EEGdenoise = pop_loadbv('.','Russell_test_2015-08-05_EEGoutsideMRIroom.vhdr') ; 
EEGdenoise = pop_chanedit(EEGdenoise,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
EEGdenoise = pop_resample(EEGdenoise,500)  ;


[es,ec] = MUK_algorithm(eeg2.data([1:end],:),64) ;
for i=1:size(ec,2) ; subplot(3,4,i) ; topoplot(ec(:,i),eeg3.chanlocs) ; end 
clear especs ; 
for i=1:size(es,1)
    [especs(i,:),f] = spectopo(es(i,:),0,500,'plot','off') ; 
end

for i=1:size(specs,1) ; subplot(5,4,i) ; plot(specs(i,f<40)) ; end

for i=1:size(EEGdenoise.data,1)
   [dspecs(i,:),f] = spectopo(EEGdenoise.data(i,:),0,500,'plot','off') ;   
   [dspecsbcg(i,:),f] = spectopo(eeg2.data(i,:),0,500,'plot','off') ; 
end

for i=1:20 ; 
    figure,
    plot(dspecs(i,f<40)) ; hold on ; plot(dspecsbcg(i,f<40),'r') ; 
    plot(specs(i,f<40),'k') ;
end

lowspecs = especs(:,f<40) ; 
lbcg = dspecs(:,f<40) ; 
c = corr(lbcg',lowspecs') ; 

[sv,si] = sort(mean(c,1),'descend') ; 
for i=1:length(si) ; subplot(4,10,i) ; plot(squeeze(especs(si(i),f<40))) ; end

for i=1:size(especs,1) ; 
    figure,
    plot(mean(dspecs,1)-15) ; hold on  ; plot(especs(i,:),'r') ; vline(20) ; 
    plot(mean(dspecsbcg,1),'k') ; 
end
    
clear espochs
for i=1:length(goodpeakinds)
    espochs(i,:,:) = es(:,goodpeakinds(i)-wprev:goodpeakinds(i)+wnext) ; 
end

stdespochs = squeeze(std(mean(espochs,1),0,3)) ; 
[sv,si] = sort(stdespochs,'descend') ; 


%%% get the FFT and subtract the mean
m1 = squeeze(mean(scalpepochs(:,1,:),1)) ; 
f1 = fft(m1) ; 
for i=1:size(scalpepochs,1)
   triali = fft(squeeze(scalpepochs(i,1,:))) ;  
   trialinv = real(ifft(triali-f1)) ; 
   subbedepochs(i,:) = trialinv ; 
end











