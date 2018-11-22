
cd c:/shared/audio/s2 ; ls 
clear all ; close all 
sets=dir('*sound*vhdr') ;
for s=1:length(sets)
   EEG = pop_loadbv('.',sets(s).name) ;  
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;

   EEG = pop_resample(EEG,500) ; 
   if s==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG) ; end    
end
merged.data = eegfiltfft(merged.data,merged.srate,1,120) ; 

merged = pop_runica(merged,'runica') ; 

triggers = {'S 12','S 13','S 14'} ;

for t=1:length(triggers)
    epochs_t = pop_epoch(merged,{triggers{t}},[-1,2.5]) ;
    for i=1:size(epochs_t.icaact,1)
        [ersp(t,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochs_t.icaact(i,:,:)),...
            epochs_t.pnts,[epochs_t.xmin,epochs_t.xmax],epochs_t.srate,0,...
            'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'baseline',0,...
            'winsize',150) ; 

    end
end

figure,for i=1:64 ; subplot(5,13,i) ; topoplot(merged.icawinv(:,i),EEG.chanlocs) ; title(i) ;  end
figure,for i=1:64 ; subplot(5,13,i)  ; imagesc(squeeze(ersp(1,i,:,:)),[-3,3]) ; title(i) ; end
figure,for i=1:64 ; subplot(5,13,i)  ; imagesc(squeeze(ersp(1,i,:,:))-squeeze(ersp(2,i,:,:)),[-3,3]) ; title(i) ; end

%%% single trials
clear sersp
for t=1:length(triggers)
    epochs_t = pop_epoch(merged,{triggers{t}},[-1,2.5]) ;
    for i=1:size(epochs_t.icaact,1)
        for k=1:size(epochs_t.icaact,3) ; 
        [sersp(t,i,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochs_t.icaact(i,:,k)),...
            epochs_t.pnts,[epochs_t.xmin,epochs_t.xmax],epochs_t.srate,0,...
            'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'baseline',0,...
            'winsize',150) ; 
        end
    end
end

clear anovas 
for a=1:20 ; disp(a) ; 
for i=1:60
    for j=1:200
        [p,anovatab,stats] = anova1(squeeze(sersp([3,1],a,:,i,j)),[],'off') ; 
        anovas(a,i,j) = anovatab{2,5} ; 
    end
end
end

figure,for i=1:20 ; subplot(4,5,i) ; imagesc(squeeze(anovas(i,:,:)),[-2,2]) ; end






