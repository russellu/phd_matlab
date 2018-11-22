clear all;
cd E:\greg_data\eeg\no_s80 ; ls
using_rest = 'true';
saving_lower_freqs = 1;
stop_before_bcg = 1;
threshold_freq = 5;

%{%
% Remove Gradient 1 (dummy volume)

% Open some sim. EEG-fMRI data and resample it.
eeg_data_loc = dir('russ_anticip_best_002.vhdr');
eeg = pop_loadbv('.', eeg_data_loc(1).name);
eeg.srate

eeg.data(32,110000:125500) = 0;

%EEG = grad_rem(EEG, 0.84, 4, 40);
tr = 0.84;
num_slices = 40;
multiband_factor = 4;
newsrate = 5000 ; 
higheeg = eeg ; 
f = fft(higheeg.data,[],2) ;  
freqincr = 5000/(size(f,2));
freqs = freqincr:freqincr:2500 ; 
freqaxis = zeros(1,size(f,2)) ; freqaxis(1:length(freqs)) = freqs ;
freqaxis(end-length(freqs)+1:end) = fliplr(freqs) ; 
raw_lowf = f(:,freqaxis<5) ; 
f(:,freqaxis<5) = 0 ; higheeg.data = real(ifft(f,[],2)) ;

highfilteeg = eegfiltfft(higheeg.data,higheeg.srate,600,5000) ; 
sliceoffset = ceil((newsrate*tr) / (num_slices/multiband_factor)) ; %375 ; %
nepochs = round(size(higheeg.data,2)/sliceoffset) ; 
gradepochs = zeros(64,nepochs,sliceoffset) ;  
highepochs = zeros(64,nepochs,sliceoffset) ; 
epochinds = zeros(nepochs,sliceoffset) ; 
icount = 1 ; 
for i=1:sliceoffset:size(higheeg.data,2)-sliceoffset
    gradepochs(:,icount,:) = higheeg.data(:,i:i+sliceoffset-1) ; 
    highepochs(:,icount,:) = highfilteeg(:,i:i+sliceoffset-1) ; 
    epochinds(icount,:) = i:i+sliceoffset-1 ; 
    icount = icount + 1 ;   
end
figure, imagesc(squeeze(gradepochs(30,:,:)));
%epochinds(epochinds==0) = 1 ; 
gradepochs = gradepochs(:,346:7535,:);
highepochs = highepochs(:,346:7535,:);
epochinds = epochinds(346:7535,:);
figure, imagesc(squeeze(gradepochs(30,:,:)));
% find the start and end (to zero pad, and remove bad epochs) 
%{
c = corr(squeeze((highepochs(1,:,:)))') ; 
c(isnan(c)) = 0 ; 
medc = median(c) ; 
grads = find(medc>.999) ; 
zc = zeros(size(medc)) ; zc(grads) = 1 ; 
% pad the non-clustered indices ; 
startclust = find(zc(1:end/2)==0) ;% startclust(end+1) = startclust(end)+1
startind = max(startclust);
size(gradepochs);
tmp = startind:startind+length(startclust)-2;
startpad = gradepochs(:,startind:startind+length(startclust)-1,:) ; 
for i=1:length(startclust) ; gradepochs(:,length(startclust)-(i-1),:) = startpad(:,i,:) ; end 

endclust = find(zc(end/2:end)==0) + round(length(zc)/2) ; 
endind = min(endclust); 
endpad = gradepochs(:,endind-length(endclust)+1:endind,:) ; 
for i=1:length(endclust) ; gradepochs(:,endclust(1)+(i-1),:) = endpad(:,end-(i-1),:) ; end
%}

subepochs = zeros(size(gradepochs)) ; 
for i=1:size(gradepochs,1) ; disp(i) ; 
    zepochs = zeros(size(gradepochs,2)+500,size(gradepochs,3)) ; 
    zepochs(251:end-250,:) = squeeze(gradepochs(i,:,:)) ; 
    zepochs(1:250,:) = zepochs(251:500,:) ; zepochs(end-250:end,:) = zepochs(end-500:end-250,:) ; 

    wsizes = [75,150,250] ; stds = wsizes/2 ;  
    filtepochs = zeros(size(zepochs)) ; 
    for w=1:length(wsizes)
    %filt = kaiser(wsizes(w)) ; filt(round(length(filt)/2)) = 0 ; filt = filt./sum(filt) ; 
    filt = fspecial('gaussian',[wsizes(w),1],stds(w)) ; filt(round(length(filt)/2)) = 0 ; filt = filt./sum(filt) ; 
    %filt = fspecial('average',[wsizes(w),1]) ; 
    filtepochs = filtepochs + imfilter(zepochs,filt)/length(wsizes) ; 
    end

    subepochs(i,:,:) = squeeze(gradepochs(i,:,:)) - filtepochs(251:end-250,:) ; 
end
neweeg = zeros(size(higheeg.data)) ; 
for i=1:size(epochinds,1)-1
   neweeg(:,epochinds(i,:)) = subepochs(:,i,:) ;     
end
for i=1:length(zc) ; if zc(i)==0 ; neweeg(:,epochinds(i,:)) = 0 ; end ; end
newhigh = higheeg ; 
newhigh.data = neweeg ; 
newhigh = pop_resample(newhigh,1024) ;