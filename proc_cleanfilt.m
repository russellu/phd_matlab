clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 
flips = [1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,0,1];
for sb=1:length(subs); disp(sb) 
    cd(['c:/shared/allres/',subs{sb}]) ; ls 
    eegs = dir('resamp*set');
    for e=1:length(eegs)
       if e==1 ; eeg = pop_loadset(eegs(e).name); else eeg = pop_mergeset(pop_loadset(eegs(e).name),eeg); end  
    end
      
    cellweights = load('cellweights'); 
    weights = cellweights.cellweights{1}; 
    sphere = cellweights.cellweights{2}; 
    winv = pinv(weights*sphere)
    subplot(4,8,sb) ; bar(zscore(mean(abs(winv(:,1:15)),2))); 
    
    bades = find(zscore(mean(abs(winv(:,1:15)),2)) > 2 | zscore(mean(abs(winv(:,1:15)),2)) < -2); 
    title([subs{sb},mat2str(bades)])
    
    cd E:\clean_allres
    mkdir(subs{sb}); 
    cd(subs{sb}); 
    save('bades','bades'); 
    flip = flips(sb); 
    save('flip','flip'); 
    pop_saveset(eeg,'eeg.set'); 
    %figure,for i=1:64 ;subplot(5,13,i)  ;topoplot(winv(:,i),eeg.chanlocs) ; title(i); end ; suptitle(subs{sb}); 
    
    
    
    %{
    filt = eegfiltfft(eeg.data,eeg.srate,1,128); 
    [weights,sphere] = runica(filt(:,1:4:end),'maxsteps',128); 
    cellweights = {weights,sphere};
    save('cellweights','cellweights') ;
    winv = pinv(weights*sphere); 
    neweeg = eeg;
    neweeg.data = weights*sphere*eeg.data; 
    epica = pop_epoch(neweeg,{'S 11','S 15','S 13'},[-.85,2.5]); 
    clear ersp; 
     for i=1:64
            [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(i,:,:)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',0,'verbose','off','timesout',100) ; 
    end
    figure,for i=1:64 ; subplot(8,16,i*2) ; imagesc(squeeze(ersp(i,:,:)),[-3,3]) ; title(i); axis xy ; colormap jet; end
    for i=1:64 ; subplot(8,16,i*2-1) ; topoplot(winv(:,i),eeg.chanlocs); end
    suptitle(subs{sb}); 
    
    %}
    %   [spectra,freq] = pwelch(x,window,overlap,Nfft,Fs,
                         %   range,plot_type,detrend,sloppy)
                         %{
   eps = pop_epoch(eeg,{'S 11','S 13','S 15'},[-1,2.5]); 

    for ch=1:64
        [base_specs,base_freqs] = pwelch(squeeze(eps.data(ch,eps.times<0,:)), 100, 50,100,eeg.srate); 
        [task_specs,task_freqs] = pwelch(squeeze(eps.data(ch,eps.times>0 & eps.times < 1000,:)), 100, 50,100,eeg.srate); 
        chanpower(ch,:) = squeeze(mean(log(task_specs) - log(base_specs),2)); 
    end
    
    allchans(sb,:,:) = chanpower; 
    %}
end

%{
subplot(1,2,1);
topoplot(squeeze(mean(mean(allchans(:,:,6),1),3)),eeg.chanlocs,'maplimits',[-0.3,0.3])
subplot(1,2,2); 
topoplot(squeeze(mean(mean(allchans(:,:,20:35),1),3)),eeg.chanlocs,'maplimits',[-0.1,0.1])
%}



