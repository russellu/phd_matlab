cd c:/shared/allres/marc ; ls 
resamps = dir('resamp_*set') ; 
for resamp=1:length(resamps)
    ecurrent = pop_loadset(resamps(resamp).name)  ;
   if resamp==1
      EEG = ecurrent ; 
   else 
      EEG = pop_mergeset(EEG,ecurrent) ; 
   end
end

epochs = pop_epoch(EEG,{'S 11','S 12','S 13','S 14','S 15','S 16'},[0,2]) ;
resdata = reshape(epochs.data,[64,size(epochs.data,2)*size(epochs.data,3)]) ; 
epochs.data = resdata ; 

epochs.data = eegfiltfft(epochs.data,epochs.srate,40,100) ; 

epochs = pop_runica(epochs,'runica') ; 


for i=1:64 ; subplot(5,13,i) ; topoplot(squeeze(epochs.icawinv(:,i)),EEG.chanlocs) ; end
copy = EEG ; resmerged = epochs ; 
copy.icaact = icaact(copy.data,resmerged.icawinv,mean(copy.data,1));
copy.icawinv = resmerged.icawinv ; copy.icasphere = resmerged.icasphere ; copy.icaweights = resmerged.icaweights ; copy.icachansind = resmerged.icachansind ; 
eps = pop_epoch(copy,{'S 14'},[-.85,2.85]) ; 

clear ersp
for i=1:64
    [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(eps.icaact(i,:,:)),size(eps.data,2),[eps.xmin,eps.xmax],eps.srate,0,'plotersp','off','plotitc','off',...
                                                    'baseline',0,'freqs',[1,150],'timesout',100,'winsize',round(eps.srate/4)) ;     
end
    
for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(i,:,:)),[-3,3]) ; end



