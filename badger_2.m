cd('C:\shared\badger\EEG-fMRI 2015-10-01\genevieve') ; 
clear all ; close all ; 
mongs = dir('*vhdr') ; 
for i=1:length(mongs) ; 
   EEG = pop_loadbv('.',mongs(i).name) ;  
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
   EEG = denoise_grad(EEG) ; 
   EEG = pop_resample(EEG,512) ; 
   EEG = denoise_bcg(EEG) ;    
   eegs{i} = EEG ; 
   if i==1 ; merged = EEG  ;else merged = pop_mergeset(EEG,merged) ; end     
end
m2 = merged ; 
filt = merged ; filt.data = eegfiltfft(filt.data,filt.srate,45,75) ; 
ica = pop_runica(filt,'runica') ; 
figure,for i=1:64 ; subplot(5,13,i) ; topoplot(squeeze(ica.icawinv(:,i)),ica.chanlocs) ; title(i) ; end

m3 = m2 ; m3.icaact = icaact(m2.data,ica.icaweights*ica.icasphere,[]) ; 
m3.icaweights = ica.icaweights ; m3.icasphere = ica.icasphere ; m3.icawinv = ica.icawinv ; m3.icachansind = ica.icachansind ; 

trigs = 1:16 ;  clear ersp ; 
for t=1:length(trigs)
    if t<10 ; trigt = ['S  ',num2str(trigs(t))] ; 
    else trigt = ['S ',num2str(trigs(t))] ; 
    end
    ep = pop_epoch(m3,{trigt},[-3,7]) ; 
    for c=1:64
        for trial=1:size(ep.data,3) ; 
            [ersp(t,c,trial,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(c,:,trial)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,'plotersp','off','plotitc','off',...
                'freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',NaN) ;
        end
    end
end

bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,200]) ;

mersp = squeeze(mean(mean(bersp,1),3)) ; 
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(times,freqs,squeeze(mersp(i,:,:)),[-8,8]) ; axis xy ; title(i) ; end ; 
comps = [24] ;
cbersp = squeeze(mean(bersp(:,comps,:,:,:),2)) ;
for i=1:16 ; figure ; for j=1:10 ; subplot(3,4,j) ; imagesc(times,freqs,squeeze(cbersp(i,j,:,:)),[-20,20]) ; axis xy ; title(['trial# ',num2str(j)]) ; end ; end

mcb = squeeze(mean(mean(cbersp(:,:,freqs>40 & freqs<70,times>0 & times<2),3),4)) ; 
mcb2 = squeeze(mean(mean(cbersp(:,:,freqs>10 & freqs<25,times>0 & times<2),3),4)) ; 
degrees = 1:23:365 ; 
errorbar(squeeze(mean(mcb,2)),squeeze(std(mcb,0,2))./sqrt(10)) ; hold on ; 
%errorbar(squeeze(mean(mcb2,2)),squeeze(std(mcb2,0,2))./sqrt(10),'r') ; 
set(gca,'XTick',1:16,'XTickLabel',degrees) ; xlabel('angle(deg)') ; title(num2str(fs)) ; 



