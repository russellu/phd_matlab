%%%% process the EEG retinotopic mapping
% quadrant starts in the bottom left, and is rotated by startangles
clear all ; close all ;
startangles = 180-[0,60,120,180,240,300] ;%-180 ;
% so stimuli are the following:
trigs{1} = {'S 11','S 12','S 13','S 14'} ; 
trigs{2} = {'S 21','S 22','S 23','S 24'} ; 
trigs{3} = {'S 31','S 32','S 33','S 34'} ; 
trigs{4} = {'S 41','S 42','S 43','S 44'} ; 
trigs{5} = {'S 51','S 52','S 53','S 54'} ; 
trigs{6} = {'S 61','S 62','S 63','S 64'} ; 

subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie'} ; 
comps= {[13,19,21,22,27,31,34],[7,15,16,17,18,21,22],[7,11,12,17,19,30],[7,23],[8,16,21],[20,27],[11,25]} ; 
for subby = 1%:7 ; 
    cd(['c:/shared/badger_eeg/',subs{subby}]) ; ls  ;
    sounds=dir('1Hz*allstim*set') ;
    for i=1:max(size(sounds)) ;  
       EEG = pop_loadset(sounds(i).name) ; 
       if i==1
          merged = EEG ; 
       else merged = pop_mergeset(EEG,merged,1) ; 
       end
    end

    clear ersp ; 
    for t=1:length(trigs)
        ep = pop_epoch(merged,trigs{t},[-2,12]) ; 
        for c=1:64 ; 
            for trial=1:size(ep.icaact,3) 
                [ersp(t,c,trial,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(c,:,trial)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                            'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'timesout',100) ; 
            end
        end
    end

    bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,size(ersp,5)]) ; 
    mbersp = squeeze(mean(bersp,3)) ; 
    for c=1:length(comps{subby}) ; figure,
        for i=1:6 ; subplot(2,3,i) ; imagesc(imfilter(squeeze(mbersp(i,comps{subby}(c),:,:)),fspecial('gaussian',3,3)),[-5,5]) ; end
    end
    
    cbersp = squeeze(mean(bersp(:,comps{subby},:,:,:),3)) ; 
    
    
end

