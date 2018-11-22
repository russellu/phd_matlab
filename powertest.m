clear all ; close all ; 

cd('C:\shared\badger\russ\russ') ; 
mongs = dir('*vhdr') ; 
for i=1:length(mongs) ; 
   EEG = pop_loadbv('.',mongs(i).name) ;  
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
   EEG = denoise_grad2(EEG) ; 
   EEG = denoise_bcg(EEG) ;    
   eegs{i} = EEG ; 
   if i==1 ; merged = EEG  ;else merged = pop_mergeset(EEG,merged) ; end     
end
mergefilt = merged ; mergefilt.data = eegfiltfft(mergefilt.data,mergefilt.srate,1,128) ; 
mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,mergefilt.srate,59.5,60.5) ; 
ica = pop_runica(mergefilt,'runica') ; 
elabs = {EEG.chanlocs.labels} ;
tp(ica) ; 
goodcs = [23,26] ; allcs = zeros(1,64) ; allcs(goodcs) = 1 ;
ica2 = pop_subcomp(ica,find(allcs==0)) ; 
comps = 1:64 ;
trigs{1} = 11:18 ; trigs{2} = 21:28 ; trigs{3} = 31:38 ; trigs{4} = 41:48 ; trigs{5} = 1:8 ; 
clear ersp 
for t1=1:length(trigs)
    for t2=1:length(trigs{t1}) ; 
        if trigs{t1}(t2) < 10
            ep = pop_epoch(ica2,{['S  ',num2str(trigs{t1}(t2))]},[-3,13]) ;
        else
            ep = pop_epoch(ica2,{['S ',num2str(trigs{t1}(t2))]},[-3,13]) ;
        end
        for c=1:length(comps) ; 
            for tr=1:size(ep.icaact,3)
                [ersp(t1,t2,c,tr,:,:),itc,powbase,times,freqs,~,~] = nonlogtimef(squeeze(ep.data(comps(c),:,tr)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                        'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'timesout',200) ; 
            end
        end
    end
end

ersp = ersp - repmat(mean(ersp(:,:,:,:,:,times<0),6),[1,1,1,1,1,200]) ; 
mersp = squeeze(mean(mean(mean(ersp,1),2),4)) ; 
fcount =1 ; 
for i=1:2:120 ; 
    subplot(6,10,fcount) ; fcount = fcount + 1 ; 
    topoplot(double(squeeze(mean(mean(mersp(:,freqs>i & freqs<i+5,times>0 & times<10),2),3))),ica.chanlocs,'maplimits',[-.1,.1]) ;
end
elabs = {ica.chanlocs.labels} ;
for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mersp(i,:,:)),[-.2,.2]) ; title(elabs{i}) ; end

