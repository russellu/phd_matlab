clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 

comps = {[5,11,14],[5,9,12],[12,21],[7,9,27],[4,5],[4,9,11,14],[3,7,12],[5,14,16],[1,14],[8,12,30],[3,6,15,19,27],[2,6,11,20],[3,21],[2,4,7],[5,7,13],[13,19],[5,7,11],[3,6,12],[5,9,12],[4,7],[4,7,13],[1,3,5],[8,12,17],[2,4,8]};


for sb=1:length(subs); disp(sb) 
    cd(['E:\clean_allres\',subs{sb}]) ; ls 
    %{
    flip = load('flip.mat'); flip = flip.flip;
    bades = load('bades.mat'); bades = bades.bades; 
    eeg = pop_loadset('eeg.set'); 
    eeg = pop_chanedit(eeg,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
    subplot(3,8,sb) ; bar(std(eeg.data,0,2)); title(subs{sb}); 
    
    newbades = bades;
    bads = find(std(eeg.data,0,2)<15);
    newbades(length(newbades)+1:length(newbades)+length(bads)) = bads; 
    save('newbades','newbades'); 
    
    if flip==1
    badvec = zeros(1,64); badvec(newbades) = 1; 
    temp_badvec = badvec(1:32) ; badvec(1:32) = badvec(33:64); badvec(33:64) = temp_badvec;
    newbades = find(badvec==1); 
    temp_eeg = eeg.data(1:32,:); eeg.data(1:32,:) = eeg.data(33:64,:); eeg.data(33:64,:) = temp_eeg; 
    end
    eeg = pop_interp(eeg,newbades,'spherical'); 
    pop_saveset(eeg,'interp_flip_eeg.set'); 
    %}
    %subplot(2,24,sb) ; topoplot(std(eeg.data,0,2),eeg.chanlocs) ; title(subs{sb}); 
    %
    %subplot(2,24,sb+24) ; topoplot(std(eeg.data,0,2),eeg.chanlocs) ; 
    
    %{
    filtdat = eegfiltfft(eeg.data,eeg.srate,1,128); 
    [weights,sphere] = runica(filtdat(:,1:5:end),'maxsteps',128);
    cleanweights{1} = weights; cleanweights{2} = sphere; 
    save('cleanweights','cleanweights'); 
    winv = pinv(weights*sphere); 

    %}
    eeg = pop_loadset('interp_flip_eeg.set'); 
    %eeg.data = eeg.data - eegfiltfft(eeg.data,eeg.srate,59,61) - eegfiltfft(eeg.data,eeg.srate,84,86); 
    cleanweights = load('cleanweights'); cleanweights = cleanweights.cleanweights; 
    weights = cleanweights{1}; sphere = cleanweights{2}; 
    winv = pinv(weights*sphere); 
    neweeg = eeg;
    acts = weights*sphere*eeg.data;     
    invacts = winv(:,comps{sb})*acts(comps{sb},:); 
    neweeg.data = invacts; 
    
    freqs = 1:2:100;
    clear log_freqtopos; 
    highfreq = eegfiltfft(invacts,eeg.srate,59,61); 
    for f=1:length(freqs)
    
    gammafilt = eegfiltfft(invacts,eeg.srate,freqs(f)-1,freqs(f)+1) - highfreq; 
    gammaeeg = eeg; gammaeeg.data = gammafilt; 
    eps = pop_epoch(gammaeeg,{'S 11','S 12','S 13','S 14','S 15','S 16'},[-1,2.5]); 
    abseps = log(abs(eps.data)); 
    baseps = abseps - repmat(mean(abseps(:,eps.times<-100,:),2),[1,size(abseps,2),1]);
    goods = find((zscore(squeeze(std(mean(baseps,1),0,2))))<2.5); 
    %subplot(3,8,sb); topoplot(squeeze(mean(mean(baseps(:,eps.times>0 & eps.times<2000,goods),2),3)),eeg.chanlocs,'maplimits',[-0.05,0.05]) ; title(subs{sb}); 
    log_freqtopos(f,:) = squeeze(mean(mean(baseps(:,eps.times>0 & eps.times<2000,goods),2),3)); 
    disp(subs{sb}); 
    %alphatopos(sb,:) = squeeze(mean(mean(baseps(:,eps.times>0 & eps.times<2000,goods),2),3)); 
    
    end
    save('log_freqtopos','log_freqtopos'); 
    %{
    epica = pop_epoch(neweeg,{'S 11','S 15','S 13'},[-.85,2.5]); 
    clear ersp; 
     for i=1:64
            [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(i,:,:)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',NaN,'verbose','off','timesout',100) ; 
     end
    figure,topoplot(squeeze(mean(mean(ersp(:,20:40,times>0 & times<2),2),3)),eeg.chanlocs,'maplimits',[-50,5])

    %figure,for i=1:64 ; subplot(8,16,i*2) ; imagesc(squeeze(ersp(i,:,:)),[-3,3]) ; title(i); axis xy ; colormap jet; end
    %for i=1:64 ; subplot(8,16,i*2-1) ; topoplot(winv(:,i),eeg.chanlocs); end
    %suptitle(subs{sb}); 
    %}
    
end

%save('alphatopos','alphatopos'); 

%{

for sb=1:length(subs); disp(sb) 
    cd(['E:\clean_allres\',subs{sb}]) ; ls 
    freqtopos = load('freqtopos');
    allfreqtopos(sb,:,:) = freqtopos.freqtopos; 
    
end

%}




