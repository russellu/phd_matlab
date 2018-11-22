cd('c:/shared/badger_eeg/valerie') ; ls  ;
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','valerie'} ; 

clear all  ; close all ;
sounds=dir('*bcg*gamma*set') ;
for i=1:max(size(sounds)) ; 
   EEG = pop_loadset(sounds(i).name) ; 
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
   alleegs{i} = EEG ; 
   if i==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG) ; end 
end

trigs = {'S  1','S  2','S  3'} ; 
clear ersp ; 
for e=1:length(alleegs)
for t=1:length(trigs)
    ep = pop_epoch(alleegs{e},{trigs{t}},[-2,7]) ; 
    for c=1:64 ; 
        for trial=1:size(ep.icaact,3)
            [ersp(t,c,(size(ep.icaact,3)*e-size(ep.icaact,3))+trial,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(c,:,trial)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                        'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'timesout',200) ; 
        end
    end
end
end
save('freqs','freqs') ; save('times','times') ; 
[s,f] = spectopo(alleegs{1}.icaact,0,256,'plot','off') ; 

bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,200]) ; 
for i=1 ; figure
    for j=1:64 
        subplot(5,13,j) ; imagesc(imfilter(squeeze(mean(bersp(i,j,:,:,:),3)),fspecial('gaussian',[3,5],3)),[-5,5]) ; title(j) ; 
    end
end
save('bersp','bersp') ; 


