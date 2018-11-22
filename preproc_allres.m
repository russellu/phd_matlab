clear all ; close all 
cd C:\shared\allres
subs = {'alex','alexandra3','audrey','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','leila','lisa','marc',...
    'marie','mathieu','maxime','menglu','mingham','olga','patricia','po','russell','suhan2','sunachakan','tah','tegan2','vincent'} ; 

badchans = {[49,54],[35,44,48],[22,42],[42],[7,17,45],[12,16,21],[17],[54,60],[54,64],[1,54],[12,35,54],[42],[22],[60,64],[54,64],[8,49],[45],[22],[13,43],...
    [16,21,54],[17],[42,45],[10,11,35,50],[34,63],[54],[17],[17,32],[51,54,60,64],[32],[8,22,42],[64]} ;

swaps = [1,1,0,0,0,0,0,1,1,1,1,0,0,1,1,1,0,0,1,1,0,0,1,0,1,0,0,1,0,0,1] ; 

for sub=1:length(subs)
    cd(['c:/shared/allres/',subs{sub}]) ; 

    
    ls
    resamps = dir('*set') ; 
    for r=1:length(resamps)
        EEG = pop_loadset(resamps(r).name) ;  
        if r==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end
    end
    if swaps(sub) ==1 
       temp = merged.data(1:32,:) ; merged.data(1:32,:) = merged.data(33:64,:) ; merged.data(33:64,:) = temp ;  
    end
    merged = pop_interp(merged,badchans{sub},'spherical') ; 
    merged = pop_resample(merged,256) ; 
    merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61) ; 
    newmerged = merged ; 
    newmerged = pop_interp(newmerged,badchans{sub},'spherical') ; 
    newmerged.data = eegfiltfft(newmerged.data,newmerged.srate,9,25) + eegfiltfft(newmerged.data,newmerged.srate,40,80)*3 ; 
    ep = pop_epoch(newmerged,{'S 11','S 13','S 14','S 15'},[-1,3]) ; 
    epica = pop_runica(ep,'runica','maxsteps',128) ;
    applmerged = merged ;
    applmerged.data = epica.icaweights*epica.icasphere*applmerged.data ; 
    
    filtweights{1} = epica.icaweights ; filtweights{2} = epica.icasphere ; 
    save('filtweights','filtweights') ; 
    
    trigs = {'S 11'} ;
    for i=1:length(trigs)
       epi = pop_epoch(applmerged,{trigs{i}},[-1,3]) ;  
       for j=1:64 ; 
          [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epi.data(j,:,:)),epi.pnts,[epi.xmin,epi.xmax],epi.srate,0,...
              'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0) ; 
       end

    end
    figure ; for j=1:64 ; subplot(5,13,j) ; imagesc(squeeze(ersp(1,j,:,:)),[-5,5]) ; end ; suptitle(subs{sub}) ; 
    
    
    
    
end



