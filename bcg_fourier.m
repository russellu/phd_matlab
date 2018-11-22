%cd('c:/Vision/Raw Files') ; ls 
cd('c:/shared/denoising_MR') ; 
clear all ; close all ; 
%EEGnoise = pop_loadbv('.','Russell_test_2015-08-05_EEGnoEPI.vhdr') ; 
%EEGnoise = pop_chanedit(EEGnoise,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
EEG = pop_loadset('merged_russell_dave.set') ; 
eeg2 = EEG ; eeg2 = pop_resample(eeg2,500) ; clear EEG

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


for i=1:size(scalpepochs,2) ; 
   meanfi = fft(squeeze(mean(scalpepochs(:,i,:),1))) ; 
   for j=1:size(scalpepochs,1) 
       fij = fft(squeeze(scalpepochs(j,i,:)))-meanfi ; 
       invs(j,i,:) = real(ifft(fij)) ; 
   end
end

alldata = eeg2.data ; 
for elec=1:64
    clear visdata2 ;
    for i=1:length(goodpeakinds)
        alldata(elec,(goodpeakinds(i))-wprev:(goodpeakinds(i))+wnext) = invs(i,elec,:) ; 
    end
end
eeg3 = eeg2 ; 
eeg3.data = alldata ; eeg3.data = eegfiltfft(eeg3.data,eeg3.srate,1,150) ; 
eeg3 = pop_runica(eeg3,'runica') ;

%eps = pop_epoch(eeg3,{'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8','S  9','S 10','S 11','S 12','S 13','S 14','S 15','S 16','S 17','S 18','S 19','S 20',},[-1,3]) ; 
%eps = pop_epoch(eeg3,{'S  1','S  4'},[-1,3]) ; 
for i=1:size(eps.icaact,1)
    [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(eps.icaact(i,:,:)),...
        eps.pnts,[eps.xmin,eps.xmax],eps.srate,0,...
        'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'baseline',0,...
        'winsize',90) ; 

end
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(i,:,:)),[-3,3]) ; end



%{
secs = 1:32*11 ; s = zeros(1,length(secs)) ; s(33:64:end) = 1 ; starts = find(s==1) ; for i=1:length(starts) ; s(starts(i):starts(i)+32) = 1 ; end 
hrf = spm_hrf(8) ; 
asl = s(1:8:end) ; 
convs = conv(asl,hrf) ; convs = convs(1:length(asl)) ; 
cd c:/shared/alex_DICOM ; dlmwrite('convs.txt',convs') ;
hrf = spm_hrf(2) ; 
bold = s(1:2:end) ; 
convs = conv(bold,hrf) ; convs = convs(1:length(bold)) ; 
cd c:/shared/alex_DICOM ; dlmwrite('convs_bold.txt',convs') ;
%}



