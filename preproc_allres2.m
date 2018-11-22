clear all ; close all 
cd C:\shared\allres
subs = {'alex','alexandra3','audrey','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','leila','lisa','marc',...
    'marie','mathieu','maxime','menglu','mingham','olga','patricia','po','russell','suhan2','sunachakan','tah','tegan2','vincent'} ; 

badchans = {[49,54],[35,44,48],[22,42],[42],[7,17,45],[12,16,21],[17],[54,60],[54,64],[1,54],[12,35,54],[42],[22],[60,64],[54,64],[8,49],[45],[22],[13,43],...
    [16,21,54],[17],[42,45],[10,11,35,50],[34,63],[54],[17],[17,32],[51,54,60,64],[32],[8,22,42],[64]} ;

swaps = [1,1,0,0,0,0,0,1,1,1,1,0,0,1,1,1,0,0,1,1,0,0,1,0,1,0,0,1,0,0,1] ; 

for sub=10:length(subs)
    cd(['c:/shared/allres/',subs{sub}]) ; 
    ls
    resamps = dir('res*set') ; 
    for r=1:length(resamps)
        EEG = pop_loadset(resamps(r).name) ;  
        if r==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end
    end
    if swaps(sub) ==1 
       temp = merged.data(1:32,:) ; merged.data(1:32,:) = merged.data(33:64,:) ; merged.data(33:64,:) = temp ;  
    end
    merged = pop_interp(merged,badchans{sub},'spherical') ; 
    merged = pop_resample(merged,256) ; 
    merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61) - eegfiltfft(merged.data,merged.srate,84,86) ; 
    filtmerged = merged ; 
    filtmerged.data = eegfiltfft(filtmerged.data,merged.srate,45,90) ;
    epmerged = pop_epoch(filtmerged,{'S 11','S 13','S 15','S 14'},[-.5,2.5]) ; 
    epica = pop_runica(epmerged,'runica') ; 
    
    
    %[weights,sphere] = runica(filtmerged.data(:,1:5:end),'maxsteps',128) ; 
    %filtmerged.data = epica.icaweights*epica.icasphere*merged.data ; winv = pinv(epica.icaweights*epica.icasphere) ; 
    filtmerged = ica_applyweights(merged,epica) ; 
    pop_saveset(filtmerged,'filtmerged') ; 
    %fullspec{1} = weights ; fullspec{2} = sphere ; 
    %save('fullspec','fullspec') ; 
    
    trigs = {'S 11'} ;
    for i=1:length(trigs)
       epi = pop_epoch(filtmerged,{trigs{i}},[-1,3]) ;  
       for j=1:64 ; 
          [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epi.icaact(j,:,:)),epi.pnts,[epi.xmin,epi.xmax],epi.srate,0,...
              'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0) ; 
       end
    end
    figure ; for j=1:64 ; subplot(5,13,j) ; imagesc(squeeze(ersp(1,j,:,:)),[-5,5]) ; end ; suptitle(subs{sub}) ; 
    %figure ; for i=1:64 ; subplot(5,13,i) ; topoplot(squeeze(winv(:,i)),merged.chanlocs) ; title(i) ; end 
    %}
end



