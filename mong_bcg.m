clear all ; close all; 
subs = {'MONG_01_RB','MONG_02_DP','MONG_03_CG','MONG_05_SG','MONG_06_TS'}; 
minpeaks = [220,160,165,150,160]; 

% plot the skew indices 
for sb=3%:length(subs)
cd(['C:\shared\mong_eeg\',subs{sb}]);
dats = dir('*dat'); 
for dat=1:length(dats)

EEG = pop_loadbv('.',strrep(dats(dat).name,'.dat','.vhdr')); 
EEG.data(32,:) = rand(1,length(EEG.data))*5; 

filt = eegfiltfft(EEG.data,EEG.srate,3,128); 
bdat = eegfiltfft(EEG.data,EEG.srate,1,128); 

[weights,sphere] = runica(filt,'maxsteps',128) ;
allweights(dat,:,:) = weights*sphere; 
acts = weights*sphere*bdat ; 
clear wskews ; 
wsize = 250*20 ; 
for i=1:5 ; jcount = 1 ; 
    for j=1:wsize:size(acts,2)-wsize
        wskews(i,jcount) = skewness(acts(i,j:j+wsize)) ; 
        jcount = jcount + 1 ; 
    end
end
skews = median((wskews),2) ; 
maxind = find(abs(skews(1:5))==max(abs(skews(1:5)))); 
if skews(maxind)<0 ; polarity = -1 ; else polarity = 1 ; end
template = acts(maxind,:)*polarity ; 
[pks,locs] = findpeaks(smooth(template),'MINPEAKDISTANCE',minpeaks(sb));
figure,plot(template); vline(locs); title(subs{sb}); 

eps = zeros(64,length(locs),200); 
epinds = zeros(length(locs),200); 
for i=2:length(locs)-1
    eps(:,i,:) = bdat(:,locs(i)-100:locs(i)+99); 
    epinds(i,:) = locs(i)-100:locs(i)+99;
end

%{
clear corrs;
for i=1:64 ; corrs(i,:,:) = corr(squeeze(eps(i,:,:))'); end
[sv,si] = sort(squeeze(mean(corrs,1)),2,'descend'); 
neweeg = EEG; 
for i=2:size(eps,2)-1
    meani = squeeze(mean(eps(:,si(3:30),:),2));     
    neweeg.data(:,epinds(i,:)) = EEG.data(:,epinds(i,:)) - meani; 
end
%}

meaneps = squeeze(mean(eps,2));
figure,plot(meaneps'); 
neweeg = EEG; 
for i=2:size(epinds,1)-1 ; 
    neweeg.data(:,epinds(i,:)) = EEG.data(:,epinds(i,:)) - meaneps; 
end
if dat==1; merged = neweeg; else merged = pop_mergeset(merged,neweeg); end
if dat==1; rawmerged = EEG; else rawmerged = pop_mergeset(rawmerged,EEG); end

end
end

filtnew = eegfiltfft(merged.data,EEG.srate,2,128); 
[weights,sphere] = runica(filtnew,'maxsteps',128) ; 
winv = pinv(weights*sphere); 
figure,for i=1:64; subplot(5,13,i) ; topoplot(winv(:,i),EEG.chanlocs) ; end
acts = weights*sphere*neweeg.data; 



