% figure_methods.m
% matlab script to preprocess EEG data and create methods figure 
% expects Brainvision denoised .set and .vhdr files
clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 
freqbands = {[46,60],[45,60],[50,70],[56,64],[50,64],[52,70],[48,60]} ; 
comps_in = {[32,29,24],[20,44],[15,29],[19,45],[24,21],[38,34],[31,25,30]} ; 
comps_out = {[17,8,4],[47,46],[49,46],[14,11],[37,19],[22,11],[45,46,48]} ; 
badtrials1 = {[],[],[],[17],[],[],[]} ; badtrials2 = {[],[],[],[31],[],[20],[]} ; badtrials3 = {[],[],[],[],[],[],[]} ; 
badouttrials = {[],[],[],[10],[],[],[]} ;
% do the ICA:
%{
for s=1:length(subs)
    cd(['c:/shared/badger_eeg/',subs{s}]) ;      
    gammas = dir('*gamma*Pulse*set') ; 
    gammas.name 
    for g=1:length(gammas)
        if g==1 ; merged = pop_loadset(gammas(g).name) ; else merged = pop_mergeset(merged,pop_loadset(gammas(g).name)) ; end
    end
    filtmerge = merged ; filtmerge.data = eegfiltfft(merged.data,merged.srate,freqbands{s}(1),freqbands{s}(2)) ; 
    filtmerge = pop_epoch(filtmerge,{'S  1','S  2','S  3'},[-1,6]) ; 
    icamerge = pop_runica(filtmerge,'runica','stop',1e-14) ; 
    epica = ica_applyweights(merged,icamerge) ; 
    pop_saveset(epica,'ica_inside.set') ; 
    epica = pop_epoch(epica,{'S  1'},[-1,6]) ; 
    for i=1:64 ; 
       [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.icaact(i,:,:)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                                                    'plotersp','off','plotitc','off','timesout',100,'freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0) ;  
        
    end
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(i,:,:)),[-6,6]) ; title(i) ; end ; suptitle(['inside ',subs(s)]) ; 

    cd outside ;  
    
    outside = dir('*outside*vhdr') ; 
    outside.name
    out = pop_loadbv('.',outside.name) ; 
    out = pop_resample(out,256) ; 
    filtout = out ; filtout.data = eegfiltfft(filtout.data,filtout.srate,freqbands{s}(1),freqbands{s}(2)) ; 
    %filtout = pop_epoch(filtout,{'S  2'},[-1,4]) ; 
    icaout = pop_runica(filtout,'runica','stop',1e-14) ; 
    epout = ica_applyweights(out,icaout) ; 
    pop_saveset(epout,'ica_outside.set') ; 
    epout = pop_epoch(epout,{'S  2'},[-1,4]) ;
     for i=1:64 ; 
       [outersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epout.icaact(i,:,:)),epout.pnts,[epout.xmin,epout.xmax],epout.srate,0,...
                                                    'plotersp','off','plotitc','off','timesout',100,'freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0) ;     
    end
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(outersp(i,:,:)),[-6,6]) ; title(i) ; end ; suptitle(['inside ',subs(s)]) ; suptitle(['outside ',subs(s)]) ;    
end
%}


for s=1:length(subs)
    cd(['c:/shared/badger_eeg/',subs{s}]) ;
    
    inside = pop_loadset('ica_inside.set') ; 
    epinside = pop_epoch(inside,{'S  1'},[-2,7]) ; 
    clear ersp ; 
    for i=1:length(comps_in{s}) ; 
        for j=1:32 ; 
           [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epinside.icaact(comps_in{s}(i),:,j)),epinside.pnts,[epinside.xmin,epinside.xmax],epinside.srate,0,...
                                                        'plotersp','off','plotitc','off','timesout',100,'freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN) ;  
        end
    end
    ersp = ersp - repmat(mean(ersp(:,:,:,times<0),4),[1,1,1,100]) ; 
    badts = zeros(1,32) ; badts(badtrials1{s}) = 1 ; goodts = (badts==0) ; 
    mersp = squeeze(mean(mean(ersp(:,goodts,:,:),1),2)) ; 
    %figure,for i=1:32 ; subplot(4,8,i) ; imagesc(squeeze(mean(ersp(:,i,:,:),1)),[-10,10]) ; title(i) ; end ; suptitle(subs{s}) ; 
    %figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(ersp(i,:,:,:),2)),[-8,8]) ; title(i) ; end
    %subplot(2,8,s) ; imagesc(mersp,[-8,8]) ; 
    allinmersp(s,:,:) = mersp ; 
    
    
    cd outside
    outside = pop_loadset('ica_outside.set') ; 
    epoutside = pop_epoch(outside,{'S  2'},[-2,5]) ; 
    clear ersp ; 
    for i=1:length(comps_out{s}) ; 
        for j=1:15
           [ersp(i,j,:,:),itc,powbase,timesout,freqsout,~,~] = newtimef(squeeze(epoutside.icaact(comps_out{s}(i),:,j)),epoutside.pnts,[epoutside.xmin,epoutside.xmax],epoutside.srate,0,...
                                                        'plotersp','off','plotitc','off','timesout',100,'freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN) ;  
        end
    end
    ersp = ersp - repmat(mean(ersp(:,:,:,timesout<0),4),[1,1,1,100]) ; 
    badts = zeros(1,15) ; badts(badouttrials{s}) = 1 ; goodts = find(badts==0) ; 
    mersp = squeeze(mean(mean(ersp(:,goodts,:,:),1),2)) ; 
    alloutmersp(s,:,:) = mersp ; 
    %figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(i,:,:)),[-8,8]) ; title(i) ; end
    %figure,for i=1:15 ; subplot(3,5,i) ; imagesc(imfilter(squeeze(mean(ersp(:,i,:,:),1)),fspecial('gaussian',[3,5],3)),[-10,10]) ; title(i) ; end ; suptitle(subs{s}) ; 
    %subplot(2,8,s+8) ; imagesc(mersp,[-8,8]) ; 
    
end

for i=1:7
   subplot(2,7,i) ; imagesc(times,freqs,squeeze(allinmersp(i,:,:)),[-6,6]) ; axis xy ; vline([0,5],'k') ; 
   subplot(2,7,i+7) ; imagesc(timesout,freqsout,squeeze(alloutmersp(i,:,:)),[-6,6]) ; axis xy ; vline([0,3],'k') ; 
   if i==1 ; xlabel('time(s)') ; ylabel('freq(hz)') ; end   
end


subplot(2,1,1) ; imagesc(times,freqs,squeeze(mean(allinmersp)),[-6,6]) ; axis xy ; vline([0,5],'k') ; 
subplot(2,1,2) ; imagesc(timesout,freqsout,squeeze(mean(alloutmersp)),[-6,6]) ; axis xy ; vline([0,3],'k') ;  xlabel('time(s)') ; ylabel('freq(hz)') ;
imagesc([-6,6]) ; colorbar ; 

subplot(1,2,1) ; 
shadedErrorBar([],squeeze(mean(mean(alloutmersp(:,:,timesout>0 & timesout<3),1),3)),squeeze(std(mean(alloutmersp(:,:,timesout>0 & timesout<3),3),0,1))./sqrt(7),'b') ;  hold on ;
shadedErrorBar([],squeeze(mean(mean(allinmersp(:,:,times>0 & times<3),1),3)),squeeze(std(mean(allinmersp(:,:,times>0 & times<3),3),0,1))./sqrt(7),'k') ; hline(0,'k') ; 
set(gca,'XTick',1:5:60,'XTickLabel',round(freqs(1:5:end))) ; xlabel('frequency(hz)') ; ylabel('power(db)') ; 
subplot(1,2,2) ; plot(squeeze(mean(mean(alloutmersp(:,freqsout>40 & freqsout < 70,timesout>0 & timesout<3),2),3)),squeeze(mean(mean(allinmersp(:,freqs>40 & freqs<70,times>0 & times<3),2),3)),'o') ; lsline ; 
[c,p] = corr([squeeze(mean(mean(alloutmersp(:,freqsout>40 & freqsout < 70,timesout>0 & timesout<3),2),3)),squeeze(mean(mean(allinmersp(:,freqs>40 & freqs<70,times>0 & times<3),2),3))]) ; 
title(['rho=',num2str(c(1,2)),', p=',num2str(p(1,2))]) ; xlabel('gamma outside (db)') ; ylabel('gamma inside (db)') ;  




subplot(2,4,1) ; topoplot(inside.icawinv(:,31),inside.chanlocs) ; subplot(2,4,5) ;imagesc(times,freqs,squeeze(ersp(31,:,:)),[-6,6]) ; axis xy ; vline([0,5],'k') ; 
subplot(2,4,2) ; topoplot(inside.icawinv(:,30),inside.chanlocs) ; subplot(2,4,6) ;imagesc(times,freqs,squeeze(ersp(30,:,:)),[-6,6]) ; axis xy ; vline([0,5],'k') ; 
subplot(2,4,3) ; topoplot(inside.icawinv(:,34),inside.chanlocs) ; subplot(2,4,7) ;imagesc(times,freqs,squeeze(ersp(34,:,:)),[-6,6]) ; axis xy ; vline([0,5],'k') ; 
subplot(2,4,4) ; topoplot(inside.icawinv(:,39),inside.chanlocs) ; subplot(2,4,8) ;imagesc(times,freqs,squeeze(ersp(39,:,:)),[-6,6]) ; xlabel('time(s)') ; ylabel('freq(hz)') ; axis xy ; vline([0,5],'k') ; 







%{
clear ersp
 for i=1:64 ; 
       [ersp(i,:,:),itc,powbase,timesout,freqsout,~,~] = newtimef(squeeze(epinside.icaact(i,:,:)),epinside.pnts,[epinside.xmin,epinside.xmax],epinside.srate,0,...
                                                    'plotersp','off','plotitc','off','timesout',100,'freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0) ;  
end
%}
