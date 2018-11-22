clear all  ;close all
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'} ; 
for s=8%:length(subs)
cd(['C:\shared\simdenoise\',subs{s}]) ;
EEG = pop_loadbv('.','retino_gamma_01.vhdr') ; 
EEG = pop_chanedit(EEG,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ; 
grad = remove_gradient2(EEG) ; 
%newgrad1 = final_bcg(grad) ; 

EEG = pop_loadbv('.','retino_gamma_02.vhdr') ; 
EEG = pop_chanedit(EEG,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ; 
grad2 = remove_gradient2(EEG) ; 
%newgrad2 = final_bcg(grad2) ; 
merged = pop_mergeset(grad2,grad) ; 
end


ep = pop_epoch(merged,{'S  1','S  2','S  3'},[-1,6]) ; %,'S  3'
lindat = reshape(ep.data,[size(ep.data,1),size(ep.data,2)*size(ep.data,3)]) ; 
%topoplot(zeros(1,64),merged.chanlocs,'electrodes','numbers') ; 
goodes = [9,20,10,60,46,31,45,59,15,51,7,37,19,38,8,52,16,29,57,27,43,23,64,24,44,28,58,30] ;%
goodes = 1:64 ; 
gooddat = lindat(goodes,:) ; gooddat = eegfiltfft(gooddat,merged.srate,40,50) ; 
[weights,sphere] = runica(gooddat,'maxsteps',128) ; 
winv = pinv(weights*sphere) ; acts = weights*sphere*lindat(goodes,:) ;
resacts = reshape(acts,[length(goodes),size(ep.data,2),size(ep.data,3)]) ; 
ep.data(1:length(goodes),:,:) = resacts ; 

%{
newmerged = merged ; newmerged.data = eegfiltfft(merged.data,merged.srate,1,60) ; 
newmergedica = pop_runica(newmerged,'runica','maxsteps',128) ; 
newmergedica.data = newmergedica.icaweights*newmergedica.icasphere*merged.data ;
ep = pop_epoch(newmergedica,{'S  1','S  2','S  3'},[-1,6]) ; %,'S  3'
%}
for j=1:length(goodes) ; 
   for k=1:32
      [dersp(j,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.data(j,:,k)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
          'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off') ; 
   end
end
dersp = dersp - repmat(mean(dersp(:,:,:,times<0),4),[1,1,1,200]) ; 
for i=1:length(goodes) ; subplot(5,13,i);  imagesc(squeeze(mean(dersp(i,:,:,:),2)),[-8,8]) ; end






%{
freqs=1:2:100 ;
actfreqs = zeros(size(acts,1),size(acts,2),length(freqs)) ; 
for i=1:length(freqs)
    actfreqs(:,:,i) = eegfiltfft(acts,merged.srate,i-1.5,i+1.5) ; 
end
resactfreqs = log(abs(reshape(actfreqs,[17,50,size(ep.data,2),size(ep.data,3)]))) ; 
baseactfreqs = resactfreqs - repmat(mean(resactfreqs(:,:,1:250,:),3),[1,1,size(resactfreqs,3),1]) ; 
%}







