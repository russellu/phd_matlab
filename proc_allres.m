cd C:\shared\allres
subs = {'alex','alexandra3','audrey','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','leila','lisa','marc',...
    'marie','mathieu','maxime','menglu','mingham','olga','patricia','po','russell','suhan2','sunachakan','tah','tegan2','vincent','ychele'} ; 

for sub=1:length(subs)
    cd(['c:/shared/allres/',subs{sub}]) ; 
    ls 
    resamps = dir('*set') ; 
    for r=1:length(resamps)
        EEG = pop_loadset(resamps(r).name) ;  
        if r==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end
    end
    newmerged = merged ; 
    newmerged.data = newmerged.data - eegfiltfft(newmerged.data,newmerged.srate,59.5,60.5) ; 
    newmerged.data = eegfiltfft(newmerged.data,merged.srate,40,80)*3 + eegfiltfft(newmerged.data,merged.srate,12,15) ; 
    eps = pop_epoch(newmerged,{'S 11','S 13','S 15','S 14'},[-1,2.2]) ; 
    icamerged = pop_runica(eps,'runica') ; 
    %tp(icamerged) ; 
    saveica{1} = icamerged.icaweights ; saveica{2} = icamerged.icasphere ; save('saveica','saveica') ; 
    applied = ica_applyweights(merged,icamerged) ; 
    trigs = {'S 11','S 13','S 15','S 14'} ;
    clear ersp ; 
    for i=1:length(trigs)
       epi = pop_epoch(applied,{trigs{i}},[-1,3]) ;  
       for j=1:64 ; 
          [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epi.icaact(j,:,:)),epi.pnts,[epi.xmin,epi.xmax],epi.srate,0,...
              'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0) ; 
       end
    end
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(ersp(1:3,i,:,:))),[-3,3]) ; title(i) ; end

end



