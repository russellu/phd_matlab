clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie'} ; 
comps = {[30,26],[32,23],[21,28],[47,46],[19,37],[47,38]} ;

for sub=1:length(subs)

    %%%% EEG processing:
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    prefix = 'highfreq_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ;gammas2=dir([prefix,'preproc*gamma*set']) ;
    allstims1=dir([prefix,'preproc*allstim*set']) ;allstims2=dir([prefix,'preproc*allstim*set']) ;
    movies=dir([prefix,'preproc*movie*set']) ; rests=dir([prefix,'preproc*rest*set']) ;
    setNames = {gammas1(1).name,gammas2(2).name} ;   
    
    trigs = {'S  1','S  2','S  3'} ; 
    
    for sn=1:length(setNames)
        EEG = pop_loadset(setNames{sn}) ; 
        for tr=1:length(trigs)
            ep = pop_epoch(EEG,trigs(tr),[-2,7]) ; 
            for i=1:length(comps{sub})
                for j=1:size(ep.icaact,3)
                    [ersp(i,tr,(sn-1)*16+j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(comps{sub}(i),:,j)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                                            'plotersp','off','plotitc','off','baseline',NaN,'timesout',200,'freqs',[1,100],'nfreqs',50,'winsize',64) ; 

                end
            end
        end
    end
   
    
   mersp = squeeze(mean(ersp,1)) ; 
   
   kersp = zeros(64,50,200) ; kersp(1:32,:,:) = squeeze(mersp(2,:,:,:)) ; kersp(33:64,:,:) = squeeze(mersp(3,:,:,:)) ; 
   
   mtkersp = squeeze(mean(kersp(:,:,times>0 & times<5),3)) ; 
   for i=1:50
      [kc,km] = kmeanscustom(uint8(mat2gray(mtkersp(:,i))*255),2) ;  
      highs = find(km==2) ; 
      lows = find(km==1) ; 
      % find the intersection between the indices and the k-means classes
      goodhigh = find(highs<=32) ; goodlow = find(lows>32) ; 
      class1(i) = length(goodhigh)/32 ; class2(i) = length(goodlow)/32 ; 
   end
   subplot(2,3,sub) ; 
   plot(class1) ; hline(0.5,'k');
   allclass1(sub,:) = class1 ; 
   allclass2(sub,:) = class2 ; 
   
end
