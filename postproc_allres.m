clear all ; close all 
cd C:\shared\allres
subs = {'alex','alexandra3','audrey','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','leila','lisa','marc',...
    'marie','mathieu','maxime','menglu','mingham','olga','patricia','po','russell','suhan2','sunachakan','tah','tegan2','vincent'} ; 

badchans = {[49,54],[35,44,48],[22,42],[42],[7,17,45],[12,16,21],[17],[54,60],[54,64],[1,54],[12,35,54],[42],[22],[60,64],[54,64],[8,49],[45],[22],[13,43],...
    [16,21,54],[17],[42,45],[10,11,35,50],[34,63],[54],[17],[17,32],[51,54,60,64],[32],[8,22,42],[64]} ;

fullcomps = {[11,12,6,13,30,38],[13,10,4,25,32,17,38,6,9],[2,5,7,8,9,15,16,22,26],[2,6,8,17,11],[5,7,6,8],[8,6,11,19,34,4],...
    [8,6,5,12],[10,14,16,22,28],[10,4,15,11],[9,5,11,12],[19,4,34],[5,16,3,26,17],[4,6,15,11,17,25],[6,13,22,29,54],[15,24,13],...
    [5,6,4,10,14],[8,1,5,13,14],[5,6,13,21],[1,13,10,18,22,23],[5,7,11,19],[4,10,5,6,8,9,29],[9,3,7,5,14],[7,6,3,17,25],[7,3,6,15,19],...
    [13,9,16],[20,18,17,3,5],[3,7,11,13,21],[8,9,5,18,44],[9,6,11,24],[9,6,8,13,15],[3,5,9,16]} ;
swaps = [1,1,0,0,0,0,0,1,1,1,1,0,0,1,1,1,0,0,1,1,0,0,1,0,1,0,0,1,0,0,1] ; 

for sub=1:length(subs)
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
    merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61) - eegfiltfft(merged.data,merged.srate,84,86) ; % remove bad harmonics

    fullspec = load('fullspec') ; fullspec = fullspec.fullspec ; 
    merged.data = fullspec{1}*fullspec{2}*merged.data ; winv = pinv(fullspec{1}*fullspec{2}) ; 
    
    bads = zeros(1,64) ; bads(fullcomps{sub}) = 1 ; badinds = find(bads==0) ; 
    newdata = merged.data ; newdata(badinds,:) = 0 ; 
    newchandata = winv*newdata ; 
    [s,f] = spectopo(newchandata,0,merged.srate,'plot','off') ; 
    figure,
    subplot(1,2,1) ; topoplot(s(:,46),merged.chanlocs) ;
    subplot(1,2,2) ; imagesc(s) ; suptitle(subs{sub}) ; 
    
    merged.data = newchandata ; 
    pop_saveset(merged,'newchandata') ; 

    %{
    trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'} ; 
    clear ersp
    for t=1:length(trigs)
        ep = pop_epoch(merged,{trigs{t}},[-.5,2.5]) ; 
        for i=1:length(fullcomps{sub}) ; 
              [ersp(t,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.data(fullcomps{sub}(i),:,:)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                  'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0) ; 
        end
    end
    
    figure,
    for i=1:3 ; subplot(2,3,i) ; imagesc(squeeze(mean(ersp(:,i,:,:))),[-3,3]) ; end 
    for i=1:3 ; subplot(2,3,i+3) ; topoplot(squeeze(winv(:,fullcomps{sub}(i))),merged.chanlocs) ; end ; suptitle(subs{sub}) ; 
    allersp(sub,:,:,:,:) = squeeze(ersp(:,1:3,:,:)) ; 
    save('stersp','ersp') ; 
    %}
    %{
    figure,
    for i=1:64 ; subplot(7,19,i) ; imagesc(squeeze(ersp(i,:,:)),[-5,5]) ; title(i) ; end
    for i=1:64 ; subplot(7,19,i+64) ; topoplot(winv(:,i),merged.chanlocs,'electrodes','off') ; title(i) ; end ;suptitle(subs{sub}) ; 
    %}
end



