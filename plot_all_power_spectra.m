clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 
compinds = {[16,25,14,12,46,64,80],[8,23,2,31,40,41],[18,40,9,15,48,61],[4,21,1,81,84,67],[53,50,14,74,79,80],[17,20,2,41,75,77],[31,33,31,48,50]} ; % fmri
epochcomps = {[31,32],[13,36],[13,26],[20,40],[14,18],[33,16],[33,37,40]} ; 
goodcomps = {[7,8,10,12,14,15,18,19,22,23,28,30],[3,8,12,13,18,25,26],[6,9,11,12,13,18,23,28,31,40],[7,8,20,24],[13,17,18,26,29,39],[8,17,18,23,32,35],[6,11,14,18,20,31,34]} ;
elecs = [59,45,31,46,60,9,20,18] ; 

for sub=1:length(subs)
    
    %%%% BOLD FMRI processing:
    cd(['c:/shared/badger_mri/',subs{sub},'/nii/melodic']) ; 
    mix = load('melodic_mix') ; 
    weights = load_untouch_nii('melodic_IC.nii.gz') ; 
    allweights{sub} = weights.img ; 
    meansub = load_untouch_nii('mean.nii.gz') ; 
    allmeans{sub} = meansub.img ; 
    mixinds = 1:735:735*5 ; 
    for i=1:length(mixinds) ; boldinds{i} = mixinds(i):mixinds(i)+734 ; end
    boldinds{6} = 735*5+1:735*5+1+449 ; 
    newboldinds{1} = boldinds{5} ; newboldinds{2} = boldinds{6} ; 
    boldinds = newboldinds ; 
    clear segmix
    for i=1:length(boldinds) ; smix = eegfiltfft(mix(boldinds{i},:)',1/0.693,0.01,2) ; segmix{i} = smix' ; end
    allsegmix{sub} = segmix ; 
    
    
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 

%{
    % get the gamma EEG power spectra (non-baseline corrected)
    gammas1=dir(['*gamma*Pulse*vhdr']) ; 
    setNames = {gammas1(1).name,gammas1(2).name} ;
    neweeg = pop_loadset('neweeg.set') ; 
    etrigs = {'S  1','S  2','S  3'} ; 
    colors = {'b','g','r'} ; 
    %subplot(2,4,sub) ; 
    for etr=1:length(etrigs)
        eps= pop_epoch(neweeg,{etrigs{etr}},[0,5]) ; 
        [s,f] = spectopo(eps.icaact,0,eps.srate,'plot','off') ; 
     %   plot(mean(s(epochcomps{sub},f>3 & f<100),1),colors{etr}) ; hold on ; 
        fs = f(f>3 & f<100) ; 
     %   set(gca,'XTick',1:50:length(fs),'XTickLabel',round(fs(1:50:end))) ; 
        specs(sub,etr,:) = mean(s(epochcomps{sub},:),1) ; 
        stimeps(sub,etr,:,:,:) = eps.icaact ; 
    end
  %  if sub==1 ; xlabel('frequency(hz)' ); ylabel('log power') ; end
    
     baseps= pop_epoch(neweeg,etrigs,[-2,0]) ; 
        [s,f] = spectopo(baseps.icaact,0,eps.srate,'plot','off') ; 
     %   plot(mean(s(epochcomps{sub},f>3 & f<100),1),colors{etr}) ; hold on ; 
        fs = f(f>3 & f<100) ; 
     %   set(gca,'XTick',1:50:length(fs),'XTickLabel',round(fs(1:50:end))) ; 
        restepcs(sub,:,:,:) = baseps.icaact ; 
    
   
    
    % plot the baseline power also
%}
  
     
    rests=dir('*rest*Pulse*set') ; 
    movies=dir('*movie*Pulse*set') ;   
    setNames = {movies(1).name,rests(1).name} ; 
    
    rcolors = {'m','k'} ; 
    weights = pop_loadset('allstim_broad_merged.set') ; 
    subplot(1,7,sub) ; 
    for setN = 1:2;    
        EEG = pop_loadset(setNames{setN}) ; 
        EEG = ica_applyweights(EEG,weights) ; 
        
        bads = ones(1,64) ; bads(goodcomps{sub}) = 0 ; 
        bads = find(bads==1) ; 
        EEG = pop_subcomp(EEG,bads) ; 
        
        lats = {EEG.urevent.latency} ; 
        types = {EEG.urevent.type} ; 
        r128s = find(strcmp('R128',types)) ; 
        rlats = cell2mat(lats(r128s)) ; 
        
        eegacts = EEG.data(:,rlats(1):rlats(length(rlats))) ; 
        eegsecs = size(eegacts,2)./EEG.srate ; 
        eegTRs = round(eegsecs./0.693)+1 ; 
  
        [s,f] = spectopo(eegacts,0,EEG.srate,'plot','off') ; 
        plot(mean(s(elecs,find(f<50))),rcolors{setN}) ; hold on ; 
        specs(sub,setN,:) = mean(s(elecs,:)) ; 
          fs = f(f>3 & f<100) ; 
        set(gca,'XTick',1:50:length(fs),'XTickLabel',round(fs(1:50:end))) ; xlim([1,length(fs)])
    end
        if sub==1 ; xlabel('frequency(hz)' ); ylabel('log power') ; end

end


freqz = find(f>1 & f<50) ; 
colors = {'m','k'} ;
for i=1:7 ; 
   subplot(2,4,i) ;  
   for j=1:2
        plot(squeeze(specs(i,j,freqz)),colors{j}) ; hold on ; 
        xlim([0,length(freqz)]) ; 
        set(gca,'XTick',1:30:length(freqz),'XTickLabel',round(f(freqz(1:30:length(freqz))))) ; 
   end
   if i==1 ; xlabel('frequency(hz)') ; ylabel('power(db)') ; end
   title(['sub ',num2str(i)]) ; 
end

shadedErrorBar([],squeeze(mean(specs(:,1,freqz),1)),squeeze(std(specs(:,1,freqz),0,1))./sqrt(7),'m') ; hold on ; 
shadedErrorBar([],squeeze(mean(specs(:,2,freqz),1)),squeeze(std(specs(:,2,freqz),0,1))./sqrt(7),'k') ; 
set(gca,'XTick',1:30:length(freqz),'XTickLabel',round(f(freqz(1:30:length(freqz))))) ; xlabel('frequency(hz)') ; ylabel('power(db)') ;


%{
for i=1:size(stimeps,1)
    for j=1:size(stimeps,2)
        for k=1:2
            for stim=1:32
                [ersp(i,j,k,stim,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(stimeps(i,j,epochcomps{i}(k),:,stim)),eps.pnts,[eps.xmin,eps.xmax],eps.srate,0,...
                    'plotersp','off','plotitc','off','timesout',50,'baseline',NaN,'winsize',120,'freqs',[1,120],'nfreqs',60) ; 
        
            end
        end
    end
end

for i=1:size(restepcs,1)
    for k=1:2
        for stim=1:32
            [bersp(i,k,stim,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(restepcs(i,epochcomps{i}(k),:,stim)),baseps.pnts,[baseps.xmin,baseps.xmax],baseps.srate,0,...
                'plotersp','off','plotitc','off','timesout',50,'baseline',NaN,'winsize',120,'freqs',[1,120],'nfreqs',60) ; 

        end
   end
end

colors = {'b','g','r'} ;
f = find(freqs>30 & freqs<90) ; freqz = freqs(f) ; 
for i=1:7 ;
   subplot(1,7,i) ; 
   stimi = squeeze(mean(mean(ersp(i,:,:,:,f,:),3),6)) ; 
   berspi = squeeze(mean(mean(bersp(i,:,:,f,:),2),5)) ; 
   for j=1:3 ; 
      shadedErrorBar([],squeeze(mean(stimi(j,:,:),2)),squeeze(std(stimi(j,:,:),0,2))./sqrt(32),colors{j}) ; 
      hold on ; 
   end
   shadedErrorBar([],squeeze(mean(berspi(:,:),1)),squeeze(std(berspi(:,:),0,1))./sqrt(32),'k') ; 
   if i==1 ; xlabel('frequency(hz)') ; ylabel('log power') ; end
   set(gca,'XTick',1:5:length(f),'XTickLabel',round(freqz(1:5:end))) ; 
   title(['subject',num2str(i)]) ; 
end

stimi = squeeze(mean(mean(mean(ersp(:,:,:,:,f,:),3),6),4)) ; 
berspi = squeeze(mean(mean(mean(bersp(:,:,:,f,:),2),5),3)) ; 

shadedErrorBar([],squeeze(mean(mean(stimi(:,1:2,:),1),2)),squeeze(std(mean(stimi(:,1:2,:),2),0,1))./sqrt(7),'b') ; hold on ; 
shadedErrorBar([],squeeze(mean(mean(stimi(:,3,:),1),2)),squeeze(std(mean(stimi(:,3,:),2),0,1))./sqrt(7),'r') ;
%shadedErrorBar([],squeeze(mean(berspi(:,:),1)),squeeze(std(berspi(:,:),0,1))./sqrt(7),'k') ; hold on ; 
xlabel('frequency(hz)') ; ylabel('log power') ;
set(gca,'XTick',1:5:length(f),'XTickLabel',round(freqz(1:5:end))) ; 
title(['grand avg',num2str(i)]) ; 
figure,plot(rand(2,3)) ; hold on ; plot(1,'k') ; legend({'0%SR','10%SR','100%SR','baseline'}) ; 
%}









%{
errorbar(squeeze(mean(specs(:,1,f>3 & f<100))),squeeze(std(specs(:,1,f>3 & f<100),0,1))./sqrt(7),'b') ;  set(gca,'XTick',1:50:length(fs),'XTickLabel',round(fs(1:50:end))) ; xlim([1,length(fs)]) ; hold on ; 
errorbar(squeeze(mean(specs(:,2,f>3 & f<100))),squeeze(std(specs(:,2,f>3 & f<100),0,1))./sqrt(7),'g') ;  set(gca,'XTick',1:50:length(fs),'XTickLabel',round(fs(1:50:end))) ; xlim([1,length(fs)])
errorbar(squeeze(mean(specs(:,3,f>3 & f<100))),squeeze(std(specs(:,3,f>3 & f<100),0,1))./sqrt(7),'r') ;  set(gca,'XTick',1:50:length(fs),'XTickLabel',round(fs(1:50:end))) ; xlim([1,length(fs)])
xlabel('frequency(hz)') ; ylabel('log power') ; 
%}

