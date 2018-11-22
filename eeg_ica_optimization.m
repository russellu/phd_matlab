subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'} ; 
for s=9%:length(subs) 
cd(['C:\shared\bscleans\',subs{s}]) ; 
gammas = dir('bsclean*gamma*set') ; 
for set=1:length(gammas)
    eeg = pop_loadset(gammas(set).name) ; 
    eegs{set} = eeg ; 
    if set==1 ; merged = eeg ; else merged = pop_mergeset(eeg,merged) ; end
end
filtmerged = eegfiltfft(merged.data,merged.srate,50,60) ;
[weights,sphere] = runica(filtmerged,'maxsteps',128) ; 
newmerged = merged ; newmerged.data = weights*sphere*merged.data ; 
allep = pop_epoch(newmerged,{'S  1','S  2'},[-2,7]) ; 
for i=1:size(allep.data,1)
    [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
            'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',200) ; 
end
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(i,:,:)),[-8,8]) ; title(i) ; end
allersp(s,:,:,:) = ersp ; 
end

comps = [20,42];  
freqs = 1:2:100 ; 
dat = weights*sphere*eegs{1}.data ; 
dat = dat(comps,:) ; clear freqdata
for f=1:length(freqs)
    freqdata(f,:,:) = eegfiltfft(dat,eeg.srate,freqs(f)-2,freqs(f)+2) ; 
end
absfreqs = log(squeeze(mean(abs(freqdata),2))) ; 
normfreqs = absfreqs - repmat(mean(absfreqs,2),[1,size(absfreqs,2)]) ; 
resfreqs = imresize(absfreqs,[50,2000]) ; 
imagesc(resfreqs(1:end,:)) ; axis xy ; 
plot(smooth(mat2gray(squeeze(mean(resfreqs(24:30,30:end-30),1)))),'r') ; hold on ; plot(smooth(mat2gray(squeeze(mean(resfreqs(5:7,30:end-30),1))))) ;
corrs = corr(resfreqs(:,100:end-100)') ; 
figure,imagesc(corrs,[-.25,.25])





