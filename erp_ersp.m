clear all ; close all  ;
cd c:/shared/resmerged ; 
subs=dir('*') ; subs(1:2) = [] ; 
for s=1:length(subs)
cd(['c:/shared/resmerged/',subs(s).name]) ; ls 
EEG = pop_loadset('merged.set') ;
newdat = eegfiltfft(EEG.data,EEG.srate,60,80)*3 +eegfiltfft(EEG.data,EEG.srate,8,25) - eegfiltfft(EEG.data,EEG.srate,59,61) ;
fulldat = EEG.data - eegfiltfft(EEG.data,EEG.srate,59,61) ;

preEEG = EEG ; preEEG.data = newdat ; preEEG = pop_epoch(preEEG,{'S 11','S 13','S 15'},[-.85,2.85]) ; 
eegica = pop_runica(preEEG,'runica','maxsteps',128) ; 
%[weights,sphere] = runica(newdat(:,1:4:end),'maxsteps',128) ; 
weights = eegica.icaweights ; sphere = eegica.icasphere ; 
winv = pinv(weights*sphere) ; 
EEG.data = weights*sphere*EEG.data ; 
epi = pop_epoch(EEG,{'S 11'},[-.85,2.5]) ; 
for i=1:size(epi.data,1) ; disp(i) ; 
[ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(epi.data(i,:,:),epi.pnts,[epi.xmin,epi.xmax],epi.srate,0,...
    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off') ; 
end
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(i,:,:)),[-5,5]) ; title(i) ; end ; suptitle(subs(s).name) ; 

saveica{1} = weights ; saveica{2} = sphere ; save('saveica','saveica') ; 

%{
comps = [14,17,18,28] ; zs = zeros(1,64) ; zs(comps) = 1 ; bads = find(zs==0) ; 
w = weights*sphere ; 
acts = weights*sphere*fulldat ; 
acts(bads,:) = 0 ; 
invacts = pinv(w)*acts ; 
specacts = invacts - eegfiltfft(invacts,EEG.srate,59,61) ; 
gammaspec = eegfiltfft(specacts,EEG.srate,8,25) ; 
EEG.data = gammaspec ; 
gammaeps = pop_epoch(EEG,{'S 11'},[-.8,2.8]) ; 
gpow = abs(gammaeps.data) ; 
gpow = gpow - repmat(mean(gpow(:,gammaeps.times<0,:),2),[1,size(gpow,2),1]) ; 
mgpow = squeeze(mean(mean(gpow(:,gammaeps.times>0 & gammaeps.times<2000,:),2),3)) ; 
topoplot(mgpow,EEG.chanlocs) ; 
%}
end

