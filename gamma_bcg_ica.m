clear all  ;close all
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'} ; 
for s=9%:length(subs)
cd(['C:\shared\simdenoise\',subs{s}]) ;

merged = pop_loadset('merged_gamma.set') ; 

%merged = newgrad ; %pop_mergeset(bgrad,bgrad2) ; 
fsubdat = eegfiltfft(merged.data,250,8,70) ;% + eegfiltfft(merged.data,250,12,14)/3 ; 
newmerged = merged ; newmerged.data = fsubdat ; 
gradep = pop_epoch(newmerged,{'S  1','S  2','S  3'},[-3,8]) ; 
gradica = pop_runica(gradep,'runica') ; 
%[nw,ns] = runica(fsubdat,'maxsteps',128) ; 
nw = gradica.icaweights ; ns = gradica.icasphere ; 
npinv = pinv(nw*ns) ; 
figure,for i=1:64 ; subplot(5,13,i) ; topoplot(npinv(:,i),merged.chanlocs) ; title(i) ; end
acts = nw*ns*merged.data ; 
newmerged = merged ; 
newmerged.data = acts ; 
trigs = {'S  1','S  2','S  3'} ; clear ersp ; 
for i=1:length(trigs)
   epi = pop_epoch(newmerged,{trigs{i}},[-2,6]) ;  
   for j=1:64 ; disp(['j=',num2str(j),', i=',num2str(i)]) ; 
       for k=1:16
      [ersp(i,j,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epi.data(j,:,k)),epi.pnts,[epi.xmin,epi.xmax],epi.srate,0,...
          'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off') ; 
       end
   end
end
bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,200]) ; 

figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(bersp(1,i,:,:,:),3)),[-8,8]) ; title(i) ; end
goodcs = [4,5,13,16,20,21,22,25,27,28,32] ;

bads = zeros(1,64) ; bads(goodcs) = 1 ; bads = find(bads==0)  ;
subacts = acts ; subacts(bads,:) = 0 ; 
newdata = npinv*subacts ; 

sumnew = sum(abs(diff(newdata,1,2)),2) ; 

% calculate ERSP on denoised electrodes
newmerged.data = newdata ; 
trigs = {'S  1','S  2','S  3'} ; clear ersp ; 
for i=1:length(trigs)
   epi = pop_epoch(newmerged,{trigs{i}},[-2,6]) ;  
   for j=1:64 ; disp(['j=',num2str(j),', i=',num2str(i)]) ; 
       for k=1:32
      [dersp(i,j,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epi.data(j,:,k)),epi.pnts,[epi.xmin,epi.xmax],epi.srate,0,...
          'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off') ; 
       end
   end
end
bdersp = dersp - repmat(mean(dersp(:,:,:,:,times<0),5),[1,1,1,1,200]) ; 
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(mean(bdersp(1:2,i,:,:,:),3),1)),[-4,4]) ; title(i) ; end

save('bdersp','bdersp') ; 



end