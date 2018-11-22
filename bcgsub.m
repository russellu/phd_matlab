clear all ; close all ; 
cd('c:/shared/badger_eeg/genevieve') ; 
grads = dir('gradient_*set') ; 
for g=1:length(grads) ; 
    
    EEG = pop_loadset(grads(g).name) ; 
    gpow = (sum(abs(EEG.data([1:31,33:end],:)),1)) ;
    [p,locs] = findpeaks(smooth(gpow,50),'MINPEAKDISTANCE',180) ; 

    clear epochs epochinds
    for i=2:length(locs)-2
        epochs(i,:,:) = EEG.data(:,locs(i)-EEG.srate/2:locs(i)+EEG.srate/2) ;
        epochinds(i,:) = locs(i)-EEG.srate/4:locs(i)+EEG.srate*(3/4) ; 
    end
    mepochs = squeeze(mean(epochs,1)) ; 
    
    clear xcorrs tshifts
    for e=1:64
    for m=1:size(epochs,1)
        xcorrs(e,m,:) = xcorr(mepochs(e,:),squeeze(epochs(m,e,:)),50) ; 
        tshifts(e,m) = find(squeeze(xcorrs(e,m,:))==max(squeeze(xcorrs(e,m,:))),1) ; 
    end
    end
    tshifts = tshifts - 50 ; 
    
    clear newepochs newepochinds
    for i=1:64
        for j=1:size(epochinds,1)
            newepochs(i,j,:) = EEG.data(i,round(epochinds(j,:)-tshifts(i,j))) ;     
            newepochinds(i,j,:) = round(epochinds(j,:)-tshifts(i,j)) ;
        end
    end
    
    clear dists ; 
    for i=1:64 ; 
        elecepochs = squeeze(newepochs(i,:,:))  ;
        dists(i,:,:) = pdist2(elecepochs,elecepochs) ; 
        for j=1:size(dists,2) ;
            [~,dists(i,j,:)] = sort(squeeze(dists(i,j,:)),'ascend') ; 
        end
    end
    
    % subtract the BCG with a moving average
    epochinds = round(newepochinds)  ;
    EEG2 = EEG ; 
    inds = 1:size(newepochinds,1) ; 
    for i=1:size(newepochinds,2)
        for j=1:64
            %avginds = find(inds>i-20 & inds < i+20 & inds ~= i) ; 
            avginds = dists(j,i,1:20) ; 
            avgepoch = squeeze(mean(newepochs(j,avginds,:),2)) ; 
            EEG2.data(j,newepochinds(j,i,:)) = EEG.data(j,newepochinds(j,i,:)) - avgepoch' ;    
        end
    end
    figure,plot(EEG.data(1,:)) ; hold on ; plot(EEG2.data(1,:),'k') ; 
    
    pop_saveset(EEG2,['bcg_',strrep(grads(g).name,'.set','')]) ; 
    
    %{
    EEG2.data = eegfiltfft(EEG2.data,EEG2.srate,5,30) ; 
    EEG2 = pop_chanedit(EEG2,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
    ica = pop_runica(EEG2,'runica') ; 
    tp(ica) ; suptitle(['g=',num2str(g)]) ;
    
    
    
    bbpow = eegfiltfft(EEG2.data,EEG2.srate,40,128).^2 ; 
    badts = (zscore(std(bbpow,0,1))>.5) ; badts = imdilate(badts,strel(ones(1,25))) ; 
    filtdat = EEG2.data ; filtdat(:,badts) = [] ;
    [weights,sphere] = runica(filtdat,'maxsteps',128) ; 
    winv = (weights*sphere)^-1 ; 
 
    figure,for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),EEG2.chanlocs) ; end ; suptitle(['g=',num2str(g)]) ;
    %}
    
    
   
    %}
    %EEG2.data = eegfiltfft(EEG2.data,EEG2.srate,1,128) ; 
    %{ 
    bbpow = eegfiltfft(EEG2.data,EEG2.srate,40,128).^2 ; 
    badts = (zscore(std(bbpow,0,1))>.5) ; badts = imdilate(badts,strel(ones(1,25))) ; 
    filtdat = EEG2.data ; filtdat(:,badts) = [] ;
    [weights,sphere] = runica(filtdat) ; 
    winv = (weights*sphere)^-1 ; 
    %}
    %{ica = pop_runica(EEG2,'runica') ; 
    %ica = pop_chanedit(ica,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
    %tp(ica) ; 
    
   % gpow = smooth(gpow,100) ;
   %{
    count = 1; 
    clear allxcs
    for c=1:1000:size(EEG.data,2)-10000
        chunk = gpow(c:c+1000) ; 
        xc(count,:) = xcorr(chunk,round(EEG.srate*4),'coeff') ; 
        cxc = xc(count,:) ; 
        plot(cxc(end/2-EEG.srate*1.5:end/2-EEG.srate*.5)) ;hold on ; 
        allxcs(count,:) = (cxc(end/2-EEG.srate*1.5:end/2-EEG.srate*.5)) ; 
        count = count + 1 ;     
    end
    %}
   %{
    for x=1:size(allxcs,1)
        hrs(x) = round(size(allxcs,2) - find((allxcs(x,:))==max((allxcs(x,:)))) + EEG.srate*.5) ;   
    end
    shr = (sin(-pi:.1:pi)) ; shr = imresize(shr,[1,hr]) ; 
    shr = repmat(shr,[1,30]) ; 
    %}
end

