cd('C:\shared\eeg_badger\MONG_01_RB\') ; 
clear all ; close all ; 
mongs = dir('*vhdr') ; 
for i=1:length(mongs) ; 
   EEG = pop_loadbv('.',mongs(i).name) ;  
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
   EEG = denoise_grad(EEG) ; 
   EEG = pop_resample(EEG,256) ; 
   EEG = denoise_bcg(EEG) ;    
   %temp = EEG.data(1:32,:) ; EEG.data(1:32,:) = EEG.data(33:64,:) ; EEG.data(33:64,:) = temp ; 
   eegs{i} = EEG ; 
   if i==1 ; merged = EEG  ;else merged = pop_mergeset(EEG,merged) ; end     
end
%diffz = zscore(sum(diff(merged.data,1,2).^2,2)) ; merged = pop_interp(merged,find(diffz>0),'spherical') ; 
filt = merged ; 
%filt.data = filt.data - eegfiltfft(filt.data,filt.srate,59.5,60.5) ; filt.data = filt.data - eegfiltfft(filt.data,filt.srate,84.5,85.5) ; 
filt.data = eegfiltfft(filt.data,filt.srate,1,128) ; 
ica = pop_runica(filt,'runica') ; 
figure,for i=1:64 ; subplot(5,13,i) ; topoplot(squeeze(ica.icawinv(:,i)),ica.chanlocs) ; title(i) ; end ; 

inds = 9 ; 
for i=1:length(inds)
    EEG = eegs{inds(i)} ; 
    if i==1 ; merged2 = EEG  ; else merged2 = pop_mergeset(EEG,merged2) ; end
end
merged2 = pop_interp(merged2,find(diffz>0),'spherical') ; 
m3 = merged2 ; m3.icaact = icaact(merged2.data,ica.icaweights*ica.icasphere,[]) ; 
m3.icaweights = ica.icaweights ; m3.icasphere = ica.icasphere ; m3.icawinv = ica.icawinv ; m3.icachansind = ica.icachansind ; 
trigs = 1:16 ;  clear ersp ; %trials = 23:24 ; 
for t=1:length(trigs)
    if t<10 ; trigt = ['S  ',num2str(trigs(t))] ; 
    else trigt = ['S ',num2str(trigs(t))] ; 
    end
    ep = pop_epoch(m3,{trigt},[-3,23]) ; 
    for c=1:64
        for trial=1:size(ep.data,3); %size(ep.data,3) ; 
            [ersp(t,c,trial,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(c,:,trial)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,'plotersp','off','plotitc','off',...
                'freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',NaN) ;
        end
    end
end

bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,200]) ; 
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(mean(bersp(:,i,:,:,:),1),3)),[-5,5]) ; title(i) ; end















%{
comp1 = squeeze(bersp(:,10,:,:,:)) ; comp2 = squeeze(bersp(:,17,:,:,:)) ; comp3 = squeeze(mean(bersp(:,[10,17,23],:,:,:),2)) ; 
f = find(freqs>50 & freqs<80) ; t = find(times>0.3 &times<10) ; 
errorbar(squeeze(mean(mean(mean(comp1(:,:,f,t),2),3),4)),squeeze(std(mean(mean(comp1(:,:,f,t),3),4),0,2))./sqrt(10)) ; 
errorbar(squeeze(mean(mean(mean(comp2(:,:,f,t),2),3),4)),squeeze(std(mean(mean(comp2(:,:,f,t),3),4),0,2))./sqrt(10)) ; 
errorbar(squeeze(mean(mean(mean(comp3(:,:,f,t),2),3),4)),squeeze(std(mean(mean(comp3(:,:,f,t),3),4),0,2))./sqrt(10)) ; figure,
errorbar(squeeze(mean(mean(mean(comp3(:,:,freqs>10 & freqs<25,t),2),3),4)),squeeze(std(mean(mean(comp3(:,:,freqs>10 & freqs<25,t),3),4),0,2))./sqrt(10)) ; 
mc3 = squeeze(mean(comp3,2)) ; 
meanf = squeeze(mean(mc3(:,freqs>10 & freqs<25,t),2)) ; 
errorbar(squeeze(mean(meanf,1)),squeeze(std(meanf,0,1))./4) ;
%}
%{
scomps = [10,17,23,53,57] ; allcomps = zeros(1,64) ; allcomps(scomps) = 1 ; badcomps = find(allcomps==0) ; 
m4 = pop_subcomp(m3,badcomps) ; m4 = pop_reref(m4,[]) ; 
trigs = 1:16 ;  clear ersp ; trials = 13:22 ; 
for t=1:length(trigs)
    if t<10 ; trigt = ['S  ',num2str(trigs(t))] ; 
    else trigt = ['S ',num2str(trigs(t))] ; 
    end
    ep = pop_epoch(m4,{trigt},[-3,13]) ; 
    for c=1:64
        for trial=1:length(trials) ; %size(ep.data,3) ; 
            [ersp(t,c,trial,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.data(c,:,trials(trial))),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,'plotersp','off','plotitc','off',...
                'freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',NaN) ;
        end
    end
end
bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,200]) ; 
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(mean(ersp(:,i,:,:,:),1),3)),[-50,0]) ; title(i) ; end
stimt = find(times>0 & times<10) ; 
mersp = squeeze(mean(ersp,3)) ; 
%}