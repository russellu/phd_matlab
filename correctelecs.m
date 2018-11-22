swapnames = {'alex','alexandra3','audrey','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','leila','lisa','marc','marie','mathieu','maxime','menglu','mingham','olga',...
            'patricia','po','russell','suhan2','sunachakan','tah','tegan2','vincent'};
swapvals = {1,1,0,0,0,0,0,1,1,1,1,0,0,1,1,1,0,0,1,1,0,0,1,0,1,0,0,1,0,0,1} ;
cd(['c:/shared/allres/','alex']) ;
EEG = pop_loadset('ica_notch85.set') ; clocs=  EEG.chanlocs ; 

for s=1:max(size(swapnames)) ; 
cd(['c:/shared/allres/',swapnames{s}]) ;
EEG = pop_loadset('ica_notch85.set') ; 
if s==32 ; EEG.chanlocs = clocs ; end 
epochs = pop_epoch(EEG,{'S 11'},[-1.000,3.000]) ; 
edat = epochs.data ; clear ersp
for i=1:size(edat,1) 
    [ersp(i,:,:),itc,powbase,times,freqs,econf,iconf] = newtimef(squeeze(edat(i,:,:)),size(edat,2),[epochs.xmin,epochs.xmax],epochs.srate,0,...
                                                            'plotersp','off','plotitc','off','baseline',0,'nfreqs',60,'freqs',[1,120]) ; 
end
if swapvals{s}
temp = ersp(1:32,:,:) ; 
ersp(1:32,:,:) = ersp(33:64,:,:) ; 
ersp(33:64,:,:) = temp ; 
end
t = find(times<2 & times>0) ; f = find(freqs<80 & freqs>50) ; 
allersp(s,:,:,:) = ersp ; 
end

 for i=1:31 ; 
    subplot(6,6,i) ; 
    topoplot(double(squeeze(mean(mean(allersp(i,:,f,t),3),4))),EEG.chanlocs) ; 
    title(swapnames{i}) ; 
end

