clear all  ; close all ; 
subs = {'alexandra3','suhan2','marie','maxime','russell','tegan2'} ;
comps = {[4,8,11,15,22,31,32],[],[],[],[],[]} ; 
trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'} ; 

clear ersp 
for s=1%:length(subs)
    cd(['C:\shared\allres\',subs{s}]) ; 
    ls
    EEG = pop_loadset('','ica_notch85.set') ; 
    for trig=1:length(trigs)
        ep = pop_epoch(EEG,{trigs{trig}},[-1,3.5]) ; 
        for c=1:64
            [ersp(s,trig,c,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(c,:,:)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                                     'plotersp','off','plotitc','off','winsize',64,'freqs',[1,128],'nfreqs',64,'baseline',0) ; 
        end
    end
end
for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(1,1,i,:,:)),[-4,4]) ; title(i) ; end

comps = [4,8,11,15,22,25,31,32,33,53] ; bads = zeros(1,64) ; bads(comps) = 1 ; 
eegsub = pop_subcomp(EEG,find(bads==0)) ; 
eegsub = pop_reref(eegsub,[]) ; 
clear stersp ; 
for trig=1:length(trigs)
    ep = pop_epoch(eegsub,{trigs{trig}},[-1,3.5]) ; 
    for c=1:64
        for trial=1:135
            [stersp(s,trig,c,trial,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.data(c,:,trial)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                                     'plotersp','off','plotitc','off','winsize',60,'freqs',[1,120],'nfreqs',60,'baseline',NaN) ; 
        end
    end
end
bersp = stersp - repmat(mean(stersp(:,:,:,:,:,times<0),6),[1,1,1,1,1,200]) ; 
mbase = squeeze(mean(stersp(:,:,:,:,:,times<0),6)) ; mtask = squeeze(mean(stersp(:,:,:,:,:,times>0 & times<2),6)) ; 
mdiff = mtask - mbase ; 


