clear all ; close all; 
cd E:\badger_eeg\russell
eeg2 = pop_loadbv('.','retino_movie.vhdr'); 
reseeg2 = pop_resample(eeg2,250); 

gradeeg = remove_gradient2(eeg2); 
gradeeg = pop_chanedit(gradeeg,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;

%{
[specs, freqs] = spectopo(gradeeg.data,0,gradeeg.srate,'plot','off'); 
[specs2, freqs2] = spectopo(reseeg2.data,0,gradeeg.srate,'plot','off'); 

plot(freqs,specs2(48,:),'r'); hold on ; plot(freqs,specs(48,:),'b'); xlabel('frequency(hz)'); ylabel('log \muV'); 
legend({'BEFORE subtraction','AFTER subtraction'});
%}


filt = eegfiltfft(gradeeg.data,gradeeg.srate,1,128); % high pass
[weights,sphere] = runica(filt(:,1:end),'maxsteps',128) ; % ICA
allweights = weights*sphere; % save weights
acts = weights*sphere*filt; 

bdat = eegfiltfft(gradeeg.data,gradeeg.srate,1,128); % for average subtraction (we are not removing values under 1Hz)
hdat = eegfiltfft(gradeeg.data,gradeeg.srate,3,128); % for individual BCG epoch correlation (more accurate to exclude low-freq)

template = acts(1,:)*-1 ; % get the bcg component time series selected in previous script
[pks,locs] = findpeaks(smooth(template),'MINPEAKDISTANCE',200); % get the peaks, use subject specific value for minpeaks

subplot(2,1,1); 
plot(mat2gray(gradeeg.data(32,33500:35000))); title('ECG channel'); xlabel('sample (HZ=250)'); ylabel('\muV');
subplot(2,1,2); 
plot(mat2gray(template(33500:35000))); title('ICA component #1'); xlabel('sample (HZ=250)'); ylabel('A.U.');

epochls = {[100,200-15]}; sb=1; 

elength = epochls{sb}(1) + epochls{sb}(2); % length of one bcg epoch (samples)
eprev = epochls{sb}(1); epost = epochls{sb}(2); % before and after peak number of samples
eps = zeros(64,length(locs),elength+1); % lower highpass epochs
heps = zeros(64,length(locs),elength+1); % higher highpass epochs
epinds = zeros(length(locs),elength+1); % epoch indices
for i=2:length(locs)-2 % get the BCG epochs using locs from findpeaks
    eps(:,i,:) = bdat(:,locs(i)-eprev:locs(i)+epost); 
    heps(:,i,:) = hdat(:,locs(i)-eprev:locs(i)+epost); 
    epinds(i,:) = locs(i)-eprev:locs(i)+epost;
end
meaneps = squeeze(mean(eps,2)); % mean lower highpass epochs (all channels, mean across epoch)

[bsv,bsi] = sort(sum(abs(meaneps),2),'descend'); % sort channels by epoch magnitude 
clear corrs
for i=1:size(heps,1) % create correlation matrix by correlating all epochs
    corrs(i,:,:) = corr(squeeze(heps(i,:,:))');  

end
mcorrs = squeeze(mean(corrs(:,:,:),1)); % mean across epochs, can use bsi to only average channels with more BCG
[sv,si] = sort(mcorrs,2,'descend') ; % sort according to highest correlation for each epoch
neweeg = gradeeg; % create new EEG for subtracted data
for i=2:size(epinds,1)- 2
    for j=1:64
        neweeg.data(j,epinds(i,:)) = squeeze(gradeeg.data(j,epinds(i,:))) - squeeze(mean(eps(j,si(i,5:40),:),2))' ; % subtract mean of top correlating (5-40) epochs
    end
end

pre_vis = eegfiltfft(gradeeg.data,gradeeg.srate,1,128); 
post_vis = eegfiltfft(neweeg.data,neweeg.srate,1,128); 

plot(pre_vis(46,25700:26800),'r') ; xlabel('time(sample), samplerate = 200Hz'); ylabel('\muV'); 
hold on ; plot(post_vis(46,25700:26800),'b'); legend('before subtraction','after subtraction'); xlabel('sample (HZ=250)'); ylabel('\muV'); 

[specs,freqs] = spectopo(gradeeg.data,0,gradeeg.srate,'plot','off'); 
[specs2,freqs2] = spectopo(neweeg.data,0,gradeeg.srate,'plot','off'); 

plot(freqs(1:100),specs(46,1:100),'r'); hold on ;plot(freqs(1:100),specs2(46,1:100),'b'); ylabel('log \muV'); xlabel('frequency(hz)');  legend('before subtraction','after subtraction');

figure,
[postweights,postsphere] = runica(post_vis,'maxsteps',128); 
postwinv = pinv(postweights*postsphere); 
for i=1:16 ; subplot(5,8,i) ; topoplot(postwinv(:,i),gradeeg.chanlocs); title(['component ',num2str(i)]); end
[preweights,presphere] = runica(pre_vis,'maxsteps',128); 
prewinv = pinv(preweights*presphere); 
for i=1:16 ; subplot(5,8,i+24) ; topoplot(prewinv(:,i),gradeeg.chanlocs); title(['component ',num2str(i)]); end


grade = [6*50 + 6*74 + 6*70 + 6*67 + 6*65 + 3*70 + 3*72 + 3*72 + 6*68 + 3*76 + 3*65 + 3*69 + 6*36 + 3*64 + 6*58 + 6*47 + 3*72 + 6*62 + 3*73 + 3*57 + 3*45 + 3*21 + 3*72 + 6*72 + 3*72 + 3*71 + 3*59]; 
normgrade = 6+6+6+6+6+3+3+3+6+3+3+3+6+3+6+6+3+6+3+3+3+3+3+6+3+3+3; 



