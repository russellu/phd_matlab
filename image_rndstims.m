clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','valerie'} ;
comps = {[32,6,48],[9,18,20],[12,10,4],[24,9,16],[8,15,52],[11,21,33],[46,33,15],[17,7,21]} ;
clear tcorrs ; 
for sub=1:length(subs)
    cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls 
    bcgicas = dir('bcgica*gamma*set') ; 
    
    clear freqepochs boldepochs
    for scan=1:2 
        cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls 
        EEG = pop_loadset(bcgicas(scan).name) ; 
        EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
        allica = load('highfreqs') ; allica = allica.highfreqs ; 
        weights = allica{1} ; sphere = allica{2} ; 
        ica = EEG ; ica.data = weights*sphere*EEG.data ; 
        winv = pinv(weights*sphere) ;
        compacts = ica.data(comps{sub},:) ; 
         stimtrigs = {'S  1','S  2','S  3'} ; clear ersp alleps
        for s=1:length(stimtrigs)
            ep = pop_epoch(ica,{stimtrigs{s}},[-1.5,6.5]) ; 
            alleps(s,:,:,:) = ep.data ; 
        end
        compeps = alleps(:,comps{sub},:,:) ; 

        clear allspecs ; wsize = 125 ; wincr = 5 ; 
        for i=1:size(compeps,1)
            for j=1:size(compeps,2)
                for k=1:size(compeps,4)
                    clear specs ; 
                    M = wsize ; Wc = .1 ; kk = (1:M-1); s = sin(Wc*kk)./ kk ; c0 = [Wc,s]; A = toeplitz(c0);
                    [V,evals] = eig(A); % Only need the principal eigenvector
                    [emax,imax] = max(abs(diag(evals)));
                    w = V(:,end-4:end) ;clear g1 
                    g1(:,1:size(w,2)) = w ; 
                    halfl = round(M/2) ; 
                    chan = squeeze(compeps(i,j,:,k)) ; icount = 1 ; 
                    for ind=halfl+1:wincr:length(chan)-M 
                        datai = chan(ind-halfl:ind+M-halfl-1) ;
                        windowed = repmat(datai,[1,size(w,2)]).*w ; 
                        f = abs(fft(windowed,[],1)) ; 
                        spec = f(1:halfl,:) ; 
                        specs(:,icount,:) = spec ; 
                        icount = icount + 1 ; 
                    end
                    allspecs(i,j,k,:,:,:) = specs ; 
                end
            end
        end
        mallspecs(scan,:,:,:,:) = squeeze(mean(mean(allspecs(:,:,:,:,:,:),2),6)) ; 
    end 
    allbase(sub,:,:,:,:,:) = mallspecs ; r
    figure,imagesc(squeeze(mean(mean(mean(mallspecs)))))
end
timesms = ep.times(halfl:wincr:end) ; 
repbase = log(allbase) - repmat(mean(log(allbase(:,:,:,:,:,1:50)),6),[1,1,1,1,1,size(allbase,6)]) ; 

for i=1:3 ; 
    subplot(1,3,i) ; 
    imagesc(timesms,1:2:125,squeeze(mean(mean(mean(mean(repbase(:,:,i,:,:,:),1),2),3),4)),[-.3,.3]) ; 
    if i==1 ; xlabel('time(ms)') ; ylabel('frequency(hz)') ; end
    axis xy ; 
end
imagesc([-3,3]) ; colorbar ; 





