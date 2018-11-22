clear all  ;close all
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'} ; 
for s=3%:length(subs)
cd(['C:\shared\simdenoise\',subs{s}]) ;
EEG = pop_loadbv('.','retino_gamma_01.vhdr') ; 
EEG = pop_chanedit(EEG,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ; 
grad = remove_gradient2(EEG) ; 
newgrad1 = final_bcg(grad) ; 

EEG = pop_loadbv('.','retino_gamma_02.vhdr') ; 
EEG = pop_chanedit(EEG,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ; 
grad = remove_gradient2(EEG) ; 
newgrad2 = final_bcg(grad) ; 
merged = pop_mergeset(newgrad1,newgrad2) ; 
%pop_saveset(merged,'merged_gamma.set') ; 
end

mergefilt = merged ; mergefilt.data = merged.data ; 
mergefilt.data = eegfiltfft(merged.data,merged.srate,1,128) ;
mergefilt = pop_runica(mergefilt,'runica','maxsteps',200) ; 
mergefilt.icaact = mergefilt.icaweights*mergefilt.icasphere*merged.data ; 
mergedata = mergefilt ; mergedata.data = mergefilt.icaact ; 
trigs = {'S  1','S  2','S  3'} ; clear ersp ; 
for i=1:length(trigs)
   epi = pop_epoch(mergedata,{trigs{i}},[-2,6]) ;  
   for j=1:64 ; disp(['j=',num2str(j),', i=',num2str(i)]) ; 
       for k=1:16
      [dersp(i,j,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epi.data(j,:,k)),epi.pnts,[epi.xmin,epi.xmax],epi.srate,0,...
          'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off') ; 
       end
   end
end
bdersp = dersp - repmat(mean(dersp(:,:,:,:,times<0),5),[1,1,1,1,200]) ; 
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(bdersp(1,i,:,:,:),3)),[-6,6]) ; title(i) ; end


%{
mersp = squeeze(mean(mean(mean(bdersp(:,[13,17,31],:,:,:),1),2),3)) ; 
ts = (mean(abs(squeeze(mersp)),2)) ;
[lv,li] = sort(ts(1:20),'descend') ; 
[hv,hi] = sort(ts(21:50),'descend') ; 
newdat = zeros(size(merged.data)) ; 
for i=1:5
    freqlow = li(i) ; lowscale = freqlow^2/100^2 ; 
    freqhigh = hi(i)+20 ; highscale = freqhigh^2/100^2 ;
    newdat = newdat + eegfiltfft(merged.data,merged.srate,freqlow-1,freqlow+1)*lowscale ; 
    newdat = newdat + eegfiltfft(merged.data,merged.srate,freqhigh-3,freqhigh+3)*highscale ; 
end
[weights,sphere] = runica(newdat,'maxsteps',128) ; 
%winv = pinv(weights*sphere) ; for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),merged.chanlocs) ; end ; 
%icamerged = merged ; icamerged.data = newdat ; icamerged = pop_epoch(mergedata,trigs,[-2,6]) ;  
%icamerged = pop_runica(icamerged,'runica','maxsteps',128) ; weights = icamerged.icaweights ; sphere = icamerged.icasphere ; 

mergedata = merged ; mergedata.data = weights*sphere*merged.data ;
trigs = {'S  1','S  2','S  3'} ; clear ersp ; 
for i=1:length(trigs)
   epi = pop_epoch(mergedata,{trigs{i}},[-2,6]) ;  
   for j=1:64 ; disp(['j=',num2str(j),', i=',num2str(i)]) ; 
       for k=1:32
      [dersp(i,j,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epi.data(j,:,k)),epi.pnts,[epi.xmin,epi.xmax],epi.srate,0,...
          'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off') ; 
       end
   end
end
bdersp = dersp - repmat(mean(dersp(:,:,:,:,times<0),5),[1,1,1,1,200]) ; 
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(bdersp(1,i,:,:,:),3)),[-6,6]) ; title(i) ; end
figure,imagesc(squeeze(mean(mean(mean(bdersp(:,[16,19,32],1:end,:,:),1),2),3)),[-4,4])
%}



