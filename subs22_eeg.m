clear all ; close all ;
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','jeremie','julie','katrine','lisa','marc','marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','vincent'} ;    
stims = {'S 11','S 12','S 13','S 14','S 15','S 16'} ;
for sub=1:length(subs) ; 
    cd(['c:/shared/allres/',subs{sub}]) ;  
    EEG = pop_loadset('ica_notch85.set') ; 
    allwinvs(sub,:,:) = EEG.icawinv ; 
    
    ags = load('allgoods') ; ags = ags.allgoods ; 
    clear ersp ; 
    for stim=1:length(stims) ;
        current = pop_epoch(EEG,{stims{stim}},[-.85,2.85]) ; 
        for comp=1:size(current.data,1)
            for trial=1:size(current.data,3)
            [ersp{stim,comp,trial}(:,:),itc,powbase,times,freqs,econf,iconf] = newtimef(squeeze(current.icaact(comp,:,trial)),current.pnts,[current.xmin,current.xmax]*1000,current.srate,...
                                                                        0,'plotitc','off','plotersp','off','freqs',[1,120],'nfreqs',60,'winsize',round(current.srate/4)) ;            
            end
        end
    end    
    save('single_trial_ersp','ersp','-v7.3') ; 
end




%{
for sub=1%:length(subs) ; 
    cd(['c:/shared/allres/',subs{sub}]) ;  
    EEG = pop_loadset('ica_notch85.set') ; 
    allwinvs(sub,:,:) = EEG.icawinv ; 
    
    ags = load('allgoods') ; ags = ags.allgoods ; 
    clear ersp ; 
    for stim=1:length(stims) ;
        current = pop_epoch(EEG,{stims{stim}},[-.85,2.85]) ; 
        for comp=1:size(current.data,1)
            [ersp(stim,comp,:,:),itc,powbase,times,freqs,econf,iconf] = newtimef(squeeze(current.icaact(comp,:,ags{stim})),current.pnts,[current.xmin,current.xmax]*1000,current.srate,...
                                                                        0,'plotitc','off','plotersp','off','freqs',[1,120],'nfreqs',60,'winsize',round(current.srate/4)) ;            
        end
    end
    allersp(sub,:,:,:,:) = ersp ; 
    
end
%}
%{
icaact = allersp ; 
for i=1:size(allersp,1) ; 
    figure; 
    for j=1:size(allersp,3)
        subplot(5,13,j); 
        imagesc(squeeze(mean(allersp(i,[1,5],j,:,:))),[-3,3]) ; 
        
    end
end
%}














