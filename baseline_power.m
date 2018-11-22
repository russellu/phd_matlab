clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie'} ; 
comps = {[30,26],[32,23],[21,28],[47,46],[19,37],[47,38]} ;

for sub=1%:length(subs)

    %%%% EEG processing:
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    prefix = 'highfreq_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ;gammas2=dir([prefix,'preproc*gamma*set']) ;
    allstims1=dir([prefix,'preproc*allstim*set']) ;allstims2=dir([prefix,'preproc*allstim*set']) ;
    movies=dir([prefix,'preproc*movie*set']) ; rests=dir([prefix,'preproc*rest*set']) ;
    setNames = {gammas1(1).name,gammas2(2).name} ;   
    
    for sn=1:length(setNames)
        EEG = pop_loadset(setNames{sn}) ; 
        ep = pop_epoch(EEG,{'S  1','S  3'},[-2,7]) ; 
        for i=1:length(comps{sub})
            for j=1:size(ep.icaact,3)
                [ersp(i,(sn-1)*32+j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(comps{sub}(i),:,j)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                                        'plotersp','off','plotitc','off','baseline',NaN,'timesout',200,'freqs',[1,100],'nfreqs',50,'winsize',64) ; 

            end
        end
    end
    
    bersp = ersp - repmat(mean(ersp(:,:,:,times<0),4),[1,1,1,200]) ; 
    mbersp = squeeze(mean(bersp,1)) ; 
    mersp = squeeze(mean(ersp,1)) ;
    
    mbase = squeeze(mean(mersp(:,:,times<0 & times>-2),3)) ; 
    mtask = squeeze(mean(mersp(:,:,times>0 & times<5),3)) ; 
    
    %subplot(2,3,sub) ; plot(std(mtask,0,1)) ; hold on ; plot(std(mbase,0,1),'r') ;
    %imagesc(corr(mbase,mtask),[-1,1]) ; 
    
   % c = corr(mbase,mtask) ;
    bandlow = 8 ; bandhigh = 25 ; 
    subplot(2,3,sub) ; 
    plot(mean(mbase(:,freqs>bandlow & freqs<bandhigh),2),mean(mtask(:,freqs>bandlow & freqs<bandhigh),2),'b.') ; 
    title(['r=',num2str(corr(mean(mbase(:,freqs>bandlow & freqs<bandhigh),2),mean(mtask(:,freqs>bandlow & freqs<bandhigh),2)))]) ; 
    
    %subplot(2,3,sub) ; imagesc(c) ; colorbar ; 
    pows = (mean(mtask(:,freqs>40 & freqs<70),2)) ; 
    [kc,km] = kmeanscustom(uint8(mat2gray(pows)*255),2) ; 
    
end












