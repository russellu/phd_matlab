clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie','tegan'} ; 
freqranges = {[48,60],[45,60],[60,80],[50,70],[50,70],[50,70],[56,70],[55,70]} ; 
lowcomps = {[25,6,48],[26,8,18],[15,16,5],[12,5],[9,7],[25,10,3],[13,7],[14,12,11]} ;
highcomps = {[20,17,8],[43,25,50],[48,51],[29,14,24],[63,27],[39,30],[32,25],[51,50,49]} ; 
medcomps = {[24,4,30],[15,4,12],[6,21,5],[8,2,17],[14,7],[9,6,11],[5,9],[15,29,13,7]} ; 

%{
for s=1:length(subs)  ;
    cd(['c:/shared/badger_eeg/',subs{s},'/outside']) ;  
    outfiles = dir('*outside*vhdr') ; 
    disp(outfiles.name) ; 
    EEG = pop_loadbv('.',outfiles(1).name) ; 
    EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
    EEG = pop_resample(EEG,256) ; 
    filtlow = EEG ; filtlow.data = eegfiltfft(EEG.data,EEG.srate,1,128) ; 
    filthigh = EEG ; filthigh.data = eegfiltfft(EEG.data,EEG.srate,freqranges{s}(1),freqranges{s}(2)) ; 
    filtmed = EEG ; filtmed.data = eegfiltfft(EEG.data,EEG.srate,7,25) ; 
    icalow = pop_runica(filtlow,'runica') ;      
    icahigh = pop_runica(filthigh,'runica') ;      
    icamed = pop_runica(filtmed,'runica') ; 
    
    applow = ica_applyweights(EEG,icalow) ; 
    epochlow = pop_epoch(applow,{'S  2'},[-2,4]) ; 
    apphigh = ica_applyweights(EEG,icahigh) ; 
    epochhigh = pop_epoch(apphigh,{'S  2'},[-2,4]) ; 
    appmed = ica_applyweights(EEG,icamed) ; 
    epochmed = pop_epoch(appmed,{'S  2'},[-2,4]) ; 
    
    pop_saveset(epochlow,'epochlow') ; pop_saveset(epochhigh,'epochhigh') ; pop_saveset(epochmed,'epochmed') ; 
    
    for i=1:64 ; 
       [ersplow(i,:,:)] = newtimef(squeeze(epochlow.icaact(i,:,:)),epochlow.pnts,[epochlow.xmin,epochlow.xmax],epochlow.srate,0,...
           'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',0,'timesout',100,'winsize',64) ;  
       [ersphigh(i,:,:)] = newtimef(squeeze(epochhigh.icaact(i,:,:)),epochhigh.pnts,[epochhigh.xmin,epochhigh.xmax],epochhigh.srate,0,...
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
    cd(['c:/shared/badger_eeg/',subs{s},'/outside']) ;  
    low = pop_loadset('epochlow.set') ; 
    high = pop_loadset('epochhigh.set') ; 
    med = pop_loadset('epochmed.set') ; 
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
    
    % for the med frequency filtered
    for i=1:length(highcomps{s}) ;
        for j=1:15
           [erspmed(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(med.icaact(medcomps{s}(i),:,j)),high.pnts,[high.xmin,high.xmax],high.srate,0,...
               'plotitc','off','plotersp','off','freqs',[1 120],'nfreqs',60,'baseline',NaN,'timesout',100,'winsize',64) ;  
        end
    end
    
    
    
    bersplow = ersplow - repmat(mean(ersplow(:,:,:,times<0),4),[1,1,1,100]) ; 
    bersphigh = ersphigh - repmat(mean(ersphigh(:,:,:,times<0),4),[1,1,1,100]) ; 
    berspmed = erspmed - repmat(mean(erspmed(:,:,:,times<0),4),[1,1,1,100]) ; 

    meanlow(s,:,:,:) = squeeze(mean(bersplow(1:2,:,:,:),1)) ; 
    meanhigh(s,:,:,:) = squeeze(mean(bersphigh(1:2,:,:,:),1)) ; 
    meanmed(s,:,:,:) = squeeze(mean(berspmed(1:2,:,:,:),1)) ; 

    figure,
    subplot(2,2,1) ; imagesc(squeeze(mean(mean(bersplow,1),2)),[-10,10]) ; 
    subplot(2,2,2) ; imagesc(squeeze(mean(mean(bersphigh,1),2)),[-10,10]) ;
    subplot(2,2,3) ; imagesc(squeeze(mean(mean(berspmed,1),2)),[-10,10])
    figure ; 
    for i=1:length(lowcomps{s}) ; subplottight(3,4,i) ; topoplot(squeeze(low.icawinv(:,lowcomps{s}(i))),low.chanlocs) ; colorbar ; end
    for i=1:length(highcomps{s}) ; subplottight(3,4,i+4) ; topoplot(squeeze(high.icawinv(:,highcomps{s}(i))),low.chanlocs) ; colorbar ; end
    for i=1:length(medcomps{s}) ; subplottight(3,4,i+8) ; topoplot(squeeze(med.icawinv(:,medcomps{s}(i))),low.chanlocs) ; colorbar ; end

end

subplot(1,2,1) ;
imagesc(squeeze(mean(mean(meanlow))),[-8,8])
subplot(1,2,2) ;
imagesc(squeeze(mean(mean(meanhigh))),[-8,8])

mmlow = squeeze(mean(meanlow(:,:,:,:,:),2)) ; mmhigh = squeeze(mean(meanhigh(:,:,:,:,:),2)) ;
for i=1:60 ; for j=1:100 ; [h,p,pi,stats] = ttest(squeeze(mmhigh(:,i,j)),squeeze(mmlow(:,i,j))) ; tersp(i,j) = stats.tstat ; persp(i,j) = p ;  end ; end

% display outside narrow band denoising results
for i=1:8 ; 
    subplot(2,8,i) ; imagesc(times,freqs,squeeze(mmlow(i,:,:)),[-10,10]) ; axis xy ;  if i==1 ; xlabel('time(s)') ; ylabel('frequency(hz)') ; end ; vline([0,3],'k') ; 
    subplot(2,8,i+8) ; imagesc(times,freqs,squeeze(mmhigh(i,:,:)),[-10,10]) ; axis xy ;  vline([0,3],'k') ; 
end

subplot(1,2,1) ; imagesc(times,freqs,tersp,[-8,8]) ; axis xy ; vline([0,3],'k') ; colorbar  ; xlabel('times') ; ylabel('frequency(hz)') ; 
subplot(1,2,2) ; imagesc(times,freqs,medfilt2(persp),[0,.05]) ; axis xy ; vline([0,3],'k') ; colorbar ; xlabel('times') ; ylabel('frequency(hz)') ;

%}





