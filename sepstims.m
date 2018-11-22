clear all ; close all 
cd c:/shared/fresh_eeg/
trigs = {'S 12','S 13','S 14','S 15','S 32','S 33','S 34','S 35'} ;
pupils = dir('stimcomps*vhdr') ; 
for p=1:length(pupils)
    EEG = pop_loadbv('.',pupils(p).name) ; 
    res = pop_resample(EEG,250) ; 
    if p==1 ; merged = res ; else merged = pop_mergeset(res,merged) ; end 
end
temp = merged.data(1:32,:) ; merged.data(1:32,:) = merged.data(33:64,:) ; merged.data(33:64,:) = temp ; 
merged = pop_chanedit(merged,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 

mergefilt = merged ; rawmerged = merged ; rawmerged.data = rawmerged.data - eegfiltfft(merged.data,250,59,61) ;
mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,250,59,61) ; 
mergefilt.data = eegfiltfft(mergefilt.data,250,7,9) ...
                 + eegfiltfft(mergefilt.data,250,15,17) ...
                 + eegfiltfft(mergefilt.data,250,23,25) ...
                 + eegfiltfft(mergefilt.data,250,31,33) ; 



epfilt = pop_epoch(mergefilt,trigs,[-.5,1.5]) ; 
mergeica = pop_runica(epfilt,'runica') ; 

applied = ica_applyweights(rawmerged,mergeica) ; 
[s,f] = spectopo(applied.icaact,0,applied.srate,'plot','off') ; 
%applied.icaact = eegfiltfft(applied.icaact,applied.srate,5,40) ; 
%applied.data = eegfiltfft(applied.data,applied.srate,5,40) ; 

clear trialspecs basespecs 
for i=1:length(trigs) ; 
    epraw = pop_epoch(applied,{trigs{i}},[-1,2]) ; 
    for j=1:120 ; 
        [trialspecs(i,j,:,:),stimf] = spectopo(squeeze(epraw.icaact(:,epraw.times>500 & epraw.times<1000,j)),0,epraw.srate,'plot','off') ; 
        [basespecs(i,j,:,:),basef] = spectopo(squeeze(epraw.icaact(:,epraw.times>-500 & epraw.times<0,j)),0,epraw.srate,'plot','off') ; 
    end
    %plot(squeeze(mean(epraw.icaact(5,:,:),3)),'Color',[i/length(trigs),0,0]) ; hold on ; 

end
for e=1:30
figure
for i=1:8
    subplot(2,4,i)  ;
    mfreqs = squeeze(trialspecs(i,:,e,:)) - squeeze(repmat(mean(basespecs(i,:,e,:),2),[1,120,1,1])) ; 
    shadedErrorBar([],mean(mfreqs,1),std(mfreqs,0,1)) ;
end
end
% paired t-test across trials? have the frequency bands as the test
% subjects
% how to best test for differences across stimulus types?







%{
compep = epraw.icaact ; 
clear allspecs ; wsize = 250 ; wincr = 5 ; 
for i=1%:size(compep,1)
    for j=1:size(compep,3)
            clear specs ; 
            M = wsize ; Wc = .1 ; kk = (1:M-1); s = sin(Wc*kk)./ kk ; c0 = [Wc,s]; A = toeplitz(c0);
            [V,evals] = eig(A); % Only need the principal eigenvector
            [emax,imax] = max(abs(diag(evals)));
            w = V(:,end-4:end) ; clear g1 
            g1(:,1:size(w,2)) = w ; 
            halfl = round(M/2) ; 
            chan = squeeze(compep(5,:,j)) ; icount = 1 ; 
            for ind=halfl+1:wincr:length(chan)-M 
                datai = chan(ind-halfl:ind+M-halfl-1)' ;
                windowed = repmat(datai,[1,size(w,2)]).*w ; 
                f = abs(fft(windowed,[],1)) ; 
                spec = f(1:halfl,:) ; 
                specs(:,icount,:) = spec ; 
                icount = icount + 1 ; 
            end
            allspecs(i,j,:,:,:) = specs ; 
    end
end
basespecs = log(allspecs) - repmat(mean(log(allspecs(:,:,:,20:40,:)),4),[1,1,1,100,1]) ; 
%}

















a = rand(50,60,70) ; 
b = a(:,:,25:70) ; inds = find(b<.5) ; 




