clear all  ;close all
cd C:\shared\simdenoise\alex
ls 

EEG = pop_loadbv('.','retino_gamma_01.vhdr') ; 
EEG = pop_chanedit(EEG,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ; 
grad = remove_gradient2(EEG) ; 
%[~,bgrad,~] = fbcg(grad) ; 
bgrad = final_bcg(grad) ; 

EEG = pop_loadbv('.','retino_gamma_02.vhdr') ; 
EEG = pop_chanedit(EEG,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ; 
grad = remove_gradient2(EEG) ; 
%[~,bgrad2,~] = fbcg(grad) ; 
bgrad2 = final_bcg(grad) ; 

merged = pop_mergeset(bgrad,bgrad2) ; 
fsubdat = eegfiltfft(merged.data,250,50,62) ; % + eegfiltfft(merged.data,250,8,16)/3 ; 
newmerged = merged ; newmerged.data = fsubdat ; 
gradep = pop_epoch(newmerged,{'S  1','S  2','S  3'},[-3,8]) ; 
gradica = pop_runica(gradep,'runica') ; 
%[nw,ns] = runica(fsubdat,'maxsteps',128) ; 
nw = gradica.icaweights ; ns = gradica.icasphere ; 
npinv = pinv(nw*ns) ; 
figure,for i=1:64 ; subplot(5,13,i) ; topoplot(npinv(:,i),EEG.chanlocs) ; end
acts = nw*ns*merged.data ; 
merged.data = acts ; 
trigs = {'S  1','S  2','S  3'} ; clear ersp ; 
for i=1:length(trigs)
   epi = pop_epoch(merged,{trigs{i}},[-2,6]) ;  
   for j=1:64 ; disp(['j=',num2str(j),', i=',num2str(i)]) ; 
       for k=1:32
      [ersp(i,j,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epi.data(j,:,k)),epi.pnts,[epi.xmin,epi.xmax],epi.srate,0,...
          'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off') ; 
       end
   end
end
bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,200]) ; 

figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(bersp(1,i,:,:,:),3)),[-8,8]) ; title(i) ; end
mb = squeeze(mean(mean(bersp(1,[33,38],:,:,:),2),3)) ; 

pop_saveset(merged,'merged.set') ; 


