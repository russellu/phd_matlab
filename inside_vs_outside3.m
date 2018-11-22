clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie','tegan'} ; 
freqranges = {[48,60],[45,60],[60,80],[50,70],[50,70],[50,70],[56,70],[55,70]} ; 
lowcomps = {[15,25],[35,30,18],[25,17],[35,20,7],[16,15],[14,33],[25,13,35],[37,14,33]} ;
highcomps = {[26,25,21],[44,37],[39,36],[19,38,50],[20,37],[34,19],[54,24,50],[43,49,46]} ; 
medcomps = {[21,20,22],[44,21,19],[21,32,14],[44,50],[29,35],[24,55],[57,35,24],[48,49,46]} ; 
%{
for s=1:length(subs)  ;
    cd(['c:/shared/badger_eeg/',subs{s}]) ;  
    outfiles = dir('*gamma*Pulse*vhdr') ; 
    %disp(outfiles.name) ; 
    EEG = pop_loadbv('.',outfiles(1).name) ; EEG2 = pop_loadbv('.',outfiles(1).name) ; EEG = pop_mergeset(EEG,EEG2) ; 
    EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
    EEG = pop_resample(EEG,250) ; 
    filtlow = EEG ; filtlow.data = eegfiltfft(EEG.data,EEG.srate,1,128) ; 
    filthigh = EEG ; filthigh.data = eegfiltfft(EEG.data,EEG.srate,freqranges{s}(1),freqranges{s}(2)) ; 
    filtmed = EEG ; filtmed.data = eegfiltfft(EEG.data,EEG.srate,30,90) ; 

    
    icalow = pop_runica(filtlow,'runica') ;      
    icahigh = pop_runica(filthigh,'runica') ;      
    icamed = pop_runica(filtmed,'runica') ;      

    applow = ica_applyweights(EEG,icalow) ; 
    epochlow = pop_epoch(applow,{'S  2'},[-2,7]) ; 
    apphigh = ica_applyweights(EEG,icahigh) ; 
    epochhigh = pop_epoch(apphigh,{'S  2'},[-2,7]) ; 
    appmed = ica_applyweights(EEG,icamed) ; 
    epochmed = pop_epoch(appmed,{'S  2'},[-2,7]) ;
    
    pop_saveset(epochlow,'epochlow') ; pop_saveset(epochhigh,'epochhigh') ; pop_saveset(epochmed,'epochmed') ; 
    
    for i=1:64 ; 
       [ersplow(i,:,:),~,~,~,~,~,~] = newtimef(squeeze(epochlow.icaact(i,:,:)),epochlow.pnts,[epochlow.xmin,epochlow.xmax],epochlow.srate,0,...
           'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',0,'timesout',100,'winsize',64) ;  
       [ersphigh(i,:,:),~,~,~,~,~,~] = newtimef(squeeze(epochhigh.icaact(i,:,:)),epochhigh.pnts,[epochhigh.xmin,epochhigh.xmax],epochhigh.srate,0,...
           'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',0,'timesout',100,'winsize',64) ;  
       [erspmed(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochmed.icaact(i,:,:)),epochmed.pnts,[epochmed.xmin,epochmed.xmax],epochmed.srate,0,...
           'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',0,'timesout',100,'winsize',64) ; 
    end
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersplow(i,:,:)),[-10,10]) ; title(i) ; end ; suptitle([subs{s},', low']) ; 
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersphigh(i,:,:)),[-10,10]) ; title(i) ; end ; suptitle([subs{s},', high']) ; 
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(erspmed(i,:,:)),[-10,10]) ; title(i) ; end ; suptitle([subs{s},', med']) ; 

end
%}



for s=1:length(subs)  ;
    cd(['c:/shared/badger_eeg/',subs{s}]) ;  
    low = pop_loadset('epochlow.set') ; 
    high = pop_loadset('epochhigh.set') ; 
    clear ersplow ersphigh
    % for the low frequency filtered data
    for i=1:length(lowcomps{s}) 
        for j=1:15
           [ersplow(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(low.icaact(lowcomps{s}(i),:,j)),low.pnts,[low.xmin,low.xmax],low.srate,0,...
               'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',NaN,'timesout',100,'winsize',64) ;   
        end
    end
    % for the high frequency filtered data
    for i=1:length(highcomps{s}) ;
        for j=1:15
           [ersphigh(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(high.icaact(highcomps{s}(i),:,j)),high.pnts,[high.xmin,high.xmax],high.srate,0,...
               'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',NaN,'timesout',100,'winsize',64) ;  
        end
    end
    
    bersplow = ersplow - repmat(mean(ersplow(:,:,:,times<0),4),[1,1,1,100]) ; 
    bersphigh = ersphigh - repmat(mean(ersphigh(:,:,:,times<0),4),[1,1,1,100]) ; 
    
    meanlow(s,:,:,:) = squeeze(mean(bersplow(1:2,:,:,:),1)) ; 
    meanhigh(s,:,:,:) = squeeze(mean(bersphigh(1:2,:,:,:),1)) ; 

    figure,
    subplot(1,2,1) ; imagesc(squeeze(mean(mean(bersplow,1),2)),[-10,10])
    subplot(1,2,2) ; imagesc(squeeze(mean(mean(bersphigh,1),2)),[-10,10]) ; suptitle(subs{s}) ; 
end



mmlow = squeeze(mean(meanlow(:,:,:,:,:),2)) ; mmhigh = squeeze(mean(meanhigh(:,:,:,:,:),2)) ;
for i=1:60 ; for j=1:100 ; [h,p,pi,stats] = ttest(squeeze(mmhigh(:,i,j)),squeeze(mmlow(:,i,j))) ; tersp(i,j) = stats.tstat ; persp(i,j) = p ;  end ; end
% display outside narrow band denoising results
for i=1:8 ; 
    subplot(2,8,i) ; imagesc(times,freqs,squeeze(mmlow(i,:,:)),[-10,10]) ; axis xy ;  if i==1 ; xlabel('time(s)') ; ylabel('frequency(hz)') ; end ; vline([0,5],'k') ; 
    subplot(2,8,i+8) ; imagesc(times,freqs,squeeze(mmhigh(i,:,:)),[-10,10]) ; axis xy ;  vline([0,5],'k') ; 
end
subplot(1,2,1) ; imagesc(times,freqs,tersp,[-8,8]) ; axis xy ; vline([0,5],'k') ; colorbar  ; xlabel('times') ; ylabel('frequency(hz)') ; 
subplot(1,2,2) ; imagesc(times,freqs,medfilt2(persp),[0,.1]) ; axis xy ; vline([0,5],'k') ; colorbar ; xlabel('times') ; ylabel('frequency(hz)') ;













