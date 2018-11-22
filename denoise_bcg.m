function denoised = denoise_bcg(EEG) 
% expects EEG with a sampling rate of 500Hz

eeg2 = EEG ;
srate = eeg2.srate ; 
% get the peaks within a long sliding window(2 seconds) for the correlation
% function
clear maxinds hrepochs hrepochs2
peakdata = eeg2.data(32,:) ;  
winlength1 = eeg2.srate*2 ; maxinds = zeros(1,size(eeg2.data,2)) ; 
wprev = srate/4 ; 
wnext = srate/1.3 ; 
for w=1:winlength1:size(eeg2.data,2)-winlength1
    windowi = peakdata(w:w+winlength1) ; 
    if max(windowi)>0 % if we are not in the zero padded region
        maxi = find(windowi==max(windowi)) ; 
        maxinds(round(w+maxi)) = 1 ; 
    end
end
hpeakinds = find(maxinds==1) ; 
for i=2:length(hpeakinds)-1
    hrepochs(i,:) = peakdata(round(hpeakinds(i)-wprev):round(hpeakinds(i)+wnext)) ; 
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
figure,imagesc(geps);

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
eeg3.data = alldata ;
denoised = eeg3 ; 

end





