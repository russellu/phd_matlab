clear all ; close all ; 
%subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'} ; 
subs = {'MONG_01_RB','MONG_02_DP','MONG_03_CG','MONG_04_AB','MONG_05_SG','MONG_06_TS'} ; 
for subby=1%:length(subs)
cd(['c:/shared/mong_eeg/',subs{subby}]) ; ls 
name = dir(['*vhdr']) ; 
for nm=1%:length(name)
eeg = pop_loadbv('.',name(nm).name) ; 
eeg = pop_chanedit(eeg,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
newsrate = 5000 ; 
higheeg = eeg ; 
resorig = pop_resample(eeg,1024) ; 
f = fft(higheeg.data,[],2) ;  
freqincr = 5000/(size(f,2));
freqs = freqincr:freqincr:2500 ; 
freqaxis = zeros(1,size(f,2)) ; freqaxis(1:length(freqs)) = freqs ;
freqaxis(end-length(freqs)+1:end) = fliplr(freqs) ; 
raw_lowf = f(:,freqaxis<5) ; 
f(:,freqaxis<5) = 0 ; higheeg.data = real(ifft(f,[],2)) ;

highfilteeg = eegfiltfft(higheeg.data,higheeg.srate,600,5000) ; 
sliceoffset = 375 ; %ceil((newsrate*0.693) / 11) ; 
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
epochinds(epochinds==0) = 1 ; 
% find the start and end (to zero pad, and remove bad epochs) 

c = corr(squeeze((highepochs(1,:,:)))') ; 
c(isnan(c)) = 0 ; 
medc = median(c) ; 
grads = find(medc>.999) ; 
zc = zeros(size(medc)) ; zc(grads) = 1 ; 
% pad the non-clustered indices ; 
startclust = find(zc(1:end/2)==0) ; startclust(end+1) = startclust(end)+1 ; 
startind = max(startclust) + 1  ; 
startpad = gradepochs(:,startind:startind+length(startclust)-1,:) ; 
for i=1:length(startclust) ; gradepochs(:,length(startclust)-(i-1),:) = startpad(:,i,:) ; end 

endclust = find(zc(end/2:end)==0) + round(length(zc)/2) ; 
endind = min(endclust)-1 ; 
endpad = gradepochs(:,endind-length(endclust)+1:endind,:) ; 
for i=1:length(endclust) ; gradepochs(:,endclust(1)+(i-1),:) = endpad(:,end-(i-1),:) ; end

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
%{
figure,
pxx = pwelch(newhigh.data(:,10000:end-10000)',2048,'power') ;
pxx(:,32) = 0 ; pxx(isnan(pxx)) = 0 ; 
subplot(1,2,1) ; plot(log(pxx(1:800,1:5:end))) ; 
subplot(1,2,2) ; imagesc(log(pxx)') ; colormap jet ; suptitle('initial subtraction') ; 
figure ; 
pxx = pwelch(resorig.data(:,10000:end-10000)',2048,'power') ;
pxx(:,32) = 0 ; pxx(isnan(pxx)) = 0 ; 
subplot(1,2,1) ; plot(log(pxx(1:800,1:5:end))) ; 
subplot(1,2,2) ; imagesc(log(pxx)') ; colormap jet ; suptitle('original data') ; 
%}
bandf = 50 ; 
fh = fft(neweeg(:,:),[],2) ; fh2 = fh ; 
fh(:,freqaxis<bandf) = 0 ; 
ih = real(ifft(fh,[],2)) ; 
hepochs = zeros(size(subepochs)) ; 
for i=1:size(epochinds,1)
    hepochs(:,i,:) = ih(:,epochinds(i,:)) ; 
end

mepochs = squeeze(mean(hepochs,2)) ; clear corrs 
for i=1:64 ; corrs(i,:) = corr(squeeze(hepochs(i,:,:))',mepochs(i,:)') ; end
c = median(corrs,1) ; 
[sv,si] = sort(c,'descend') ; 
indsi = 1:length(si) ; 
subhepochs = zeros(size(hepochs)) ; 
for i=1:length(si) ; disp(i) ; 
    closeinds = find(abs(indsi-i)<10) ; 
    subhepochs(:,si(i),:) = squeeze(mean(hepochs(:,si(closeinds),:),2)) ; 
end
highsubs = zeros(size(eeg.data)) ; for i=1:size(epochinds,1) ; highsubs(:,epochinds(i,:)) = hepochs(:,i,:) - subhepochs(:,i,:) ; end ; 
highf = fft(highsubs,[],2) ; 
plot(abs(fh(13,:))) ; hold on ; plot(abs(highf(13,:)),'r') ;
highf(:,freqaxis<bandf) = fh2(:,freqaxis<bandf) ; highf(:,freqaxis<5) = raw_lowf ; 
iclean = real(ifft(highf,[],2)) ; 
eeg2 = eeg ; eeg2.data = iclean ; eeg2 = pop_resample(eeg2,1024) ; 
pxx = pwelch(eeg2.data(:,10000:end-10000)',2048,'power') ;
pxx(:,32) = 0 ; pxx(isnan(pxx)) = 0 ; 
figure,subplot(1,2,1) ; plot(log(pxx(1:800,1:5:end))) ; 
subplot(1,2,2) ; imagesc(log(pxx)') ; colormap jet ; suptitle(['final result ',subs{subby}]) 


%pop_saveset(eeg2,['newgrad_',strrep(name(nm).name,'.vhdr','.set')]) ; 

end
end