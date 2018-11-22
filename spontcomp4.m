clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 

comps = {[30,26],[32,23],[21,28],[47,46],[19,37],[47,38],[44,46]} ; % all subjects, right and left
fcomps = {[64,80,12],[41,31,40],[61,48,15],[6,84,79],[79,80,50],[75,77,41],[48,50]} ; 

for sub=1:length(subs)

    %%%% BOLD FMRI processing:
    cd(['c:/shared/badger_mri/',subs{sub},'/nii/melodic']) ; 
    mix = load('melodic_mix') ; 
    weights = load_untouch_nii('melodic_IC.nii.gz') ; 
    allweights{sub} = weights.img ; 
    meansub = load_untouch_nii('mean.nii.gz') ; 
    allmeans{sub} = meansub.img ; 
    mixinds = 1:735:735*5 ; 
    for i=1:length(mixinds) ; boldinds{i} = mixinds(i):mixinds(i)+734 ; end
    boldinds{6} = 735*5+1:735*5+1+449 ; 
    for i=1:length(boldinds) ; segmix{i} = mix(boldinds{i},:) ; end
    
    %%%% EEG processing:
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    prefix = 'highfreq_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ;gammas2=dir([prefix,'preproc*gamma*set']) ;
    allstims1=dir([prefix,'preproc*allstim*set']) ;allstims2=dir([prefix,'preproc*allstim*set']) ;
    movies=dir([prefix,'preproc*movie*set']) ; rests=dir([prefix,'preproc*rest*set']) ;
    setNames = {allstims1(1).name,allstims1(2).name,gammas1(1).name,gammas1(2).name,movies(1).name,rests(1).name} ;   

    for setN = 1:6;    
        EEG = pop_loadset(setNames{setN}) ; 
        lats = {EEG.urevent.latency} ; 
        types = {EEG.urevent.type} ; 
        r128s = find(strcmp('R128',types)) ; 
        rlats = cell2mat(lats(r128s)) ; 
        
        eegacts = EEG.icaact(:,rlats(1):rlats(length(rlats))) ; 
        eegsecs = size(eegacts,2)./EEG.srate ; 
        eegTRs = round(eegsecs./0.693) ; 
    
        clear epow
        for c=1:length(comps{sub})  
            [epow(c,:,:),itc,powbase,times,freqs,~,~] = newtimef(eegacts(comps{sub}(c),:),size(eegacts,2),[0,size(eegacts,2)./EEG.srate],EEG.srate,0,...
                'plotersp','off','plotitc','off','baseline',NaN,'timesout',eegTRs,'nfreqs',50,'winsize',250,'freqs',[1,100]) ; 
        end
        indlags = -20:20 ; 
        mepow = squeeze(mean(epow,1)) ; 
        smix = segmix{setN} ; smix = eegfiltfft(smix',1/0.693,0.03,1.5) ; smix = smix' ; 
        smix = mean(smix(1:size(mepow,2),fcomps{sub}(1:2)),2) ; 
        if setN == 6 
            restmix(sub,:) = smix ; 
            restpow(sub,:,:) = mepow ; 
        end 
        if setN == 5  
            moviemix(sub,:) = smix ; 
            moviepow(sub,:,:) = mepow ; 
        end 
        for i=1:length(indlags) % for all time lags
            for j=1:size(mepow,1) % for all freqs
                powj = smooth(mean(mepow(freqs>freqs(j)-2.5 & freqs<freqs(j)+2.5,:))) ; 
                corrs(j,i) = corr2(powj(50:end-50)',smix(50+indlags(i):end-50+indlags(i))') ;               
            end   
        end
        allcorrs(sub,setN,:,:) = corrs ; 
        %figure,
        %mepow = squeeze(mean(epow,1)) ; 
        %mepow = medfilt2(mepow - repmat(mean(mepow,2),[1,size(mepow,2)])) ; 
        %imagesc(times,freqs,mepow) ;
    end
end

%for i=1:6 ; subplot(2,3,i) ; imagesc(squeeze(allcorrs(i,:,:)),[-.5,.5]) ; axis xy ; end 
%figure,imagesc(squeeze(mean(allcorrs,1)),[-.3,.3]) ; axis xy ; 
subplot(2,2,1) ; imagesc(squeeze(mean(mean(allcorrs(:,1:2,:,:),1),2)),[-.35,.35]) ; axis xy ;
subplot(2,2,2) ; imagesc(squeeze(mean(mean(allcorrs(:,3:4,:,:),1),2)),[-.35,.35]) ; axis xy ;
subplot(2,2,3) ; imagesc(squeeze(mean(mean(allcorrs(:,5,:,:),1),2)),[-.35,.35]) ; axis xy ;
subplot(2,2,4) ; imagesc(squeeze(mean(mean(allcorrs(:,6,:,:),1),2)),[-.35,.35]) ; axis xy ;

tlags = -20:20 ; tlags = tlags.*0.693 ; 

tkern = squeeze(mean(mean(allcorrs(:,5,:,:),1),2)) ; 
%tkern = tkern.*(abs(tkern)>.01) ;
tkern = fliplr(tkern) ; 

% for the resting
clear pred
for s=1:7
    for i=21:size(restmix,2)-40
       pred(i) = sum(sum(squeeze(restpow(s,:,i-20:i+20)).*tkern)) ; 
    end
    subplot(7,2,(s-1)*2+1) ; 
    plot(mat2gray(pred(22:end)),mat2gray(restmix(s,22:end-40)),'r.') ; lsline
    title(['r=',num2str(corr2(mat2gray(pred(22:end)),mat2gray(restmix(s,22:end-40))))]) ; 
    subplot(7,2,(s-1)*2+2) ;  plot(mat2gray(pred(22:end))) ; hold on ; plot(mat2gray(restmix(s,22:end-40)),'k') ; 
    title(['r^2=',num2str(corr2(mat2gray(pred(22:end)),mat2gray(restmix(s,22:end-40))).^2)]) ; 
    allrsqr(s) = corr2(mat2gray(pred(22:end)),mat2gray(restmix(s,22:end-40))).^2 ; 
end
suptitle(['mean r^2=',num2str(mean(allrsqr))]) ; 

% for the movie
clear pred
for s=1:7
    for i=21:size(moviemix,2)-40
       pred(i) = sum(sum(squeeze(moviepow(s,:,i-20:i+20)).*tkern)) ; 
    end
    subplot(7,2,(s-1)*2+1) ; 
    plot(mat2gray(pred(22:end)),mat2gray(moviemix(s,22:end-40)),'r.') ; lsline
    title(['r=',num2str(corr2(mat2gray(pred(22:end)),mat2gray(moviemix(s,22:end-40))))]) ; 
    subplot(7,2,(s-1)*2+2) ;  plot(mat2gray(pred(22:end))) ; hold on ; plot(mat2gray(moviemix(s,22:end-40)),'k') ; 
    title(['r^2=',num2str(corr2(mat2gray(pred(22:end)),mat2gray(moviemix(s,22:end-40))).^2)]) ; 
    allrsqr(s) = corr2(mat2gray(pred(22:end)),mat2gray(moviemix(s,22:end-40))).^2 ; 
end
suptitle(['mean r^2=',num2str(mean(allrsqr))]) ; 



figure,
for i=1:7 ; subplot(2,7,i) ;
    imagesc(tlags,freqs,squeeze(allcorrs(i,6,:,:)),[-.5,.5]) ; axis xy ; vline(0,'k') ; 
    if i==1 ; xlabel('tlag(s)') ; ylabel('frequency(hz)') ; end
    title(['sub=',num2str(i)]) ; 
end
for i=1:7 ; subplot(2,7,i+7) ;
    imagesc(tlags,freqs,squeeze(allcorrs(i,5,:,:)),[-.5,.5]) ; axis xy ; vline(0,'k') ; 
    if i==1 ; xlabel('tlag(s)') ; ylabel('frequency(hz)') ; end
end

mc1 = squeeze(mean(allcorrs(:,5,:,:))) ; rmc1 = reshape(mc1,[1,numel(mc1)]) ; 
mc2 = squeeze(mean(allcorrs(:,6,:,:))) ; rmc2 = reshape(mc2,[1,numel(mc2)]) ; 
plot(rmc1,rmc2,'.') ; [c,p] = corr([rmc1;rmc2]') ; 
title(['rho = ',num2str(c(1,2)), ' p = ',num2str(p(1,2))]) ; xlabel('movie') ; ylabel('rest') ; lsline ; xlim([-0.3,0.2]) ; ylim([-0.3,0.2]) ;  




for i=1:7 ; subplot(1,7,i) ; imagesc(tlags,freqs,squeeze(allcorrs(i,6,:,:)),[-.4,.4]) ; axis xy ;xlabel('time(s)') ; ylabel('freq(hz)') ;vline(0,'k') ;  end ; 




