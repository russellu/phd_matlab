clear all ; close all; 
%{
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'} ;
comps = {[4,5,6,7,8,9,10,11,12,13,14,15,17,20,21,22,23,25,26,28,30,38],[5,7,8,9,13,17,18,19,20,25,26,28,29,37,47,58],...
         [4,5,6,7,8,9,11,13,14,15,16,17,20,21,23,25,26,27,28,35],[3,5,8,9,10,12,17,18,19,21,22,27,29,36],...
         [4,8,9,11,12,13,14,16,20,21,25,27,36],[5,6,7,8,9,10,12,13,14,16,20,26,35,37],...
         [1,6,8,9,11,13,15,17,24,25,28,33,35,39],[3,5,7,8,10,15,16,20,21,23,26,27,28,29,30,38,40],...
         [3,6,8,10,11,16,19,20,22,25,30,36,42]};
%}
%{
subs = {'MONG_01_RB','MONG_02_DP','MONG_03_CG','MONG_05_SG','MONG_06_TS'}; 
comps = {[3,5,6,7,8,10,11,14,16,17,21,22,24,25,26,30,31],[2,4,9,10,11,15,16,17,18,20,22,27,29,42,48],...
         [2,4,7,9,11,12,15,16,18,23,24,29,30,32,33,34,36,39,40,41,46],[5,6,7,8,9,10,11,15,16,17,18,19,20,24,25,37,48,49],...
         [6,10,13,15,16,18,19,21,23,25,28,41,43,52,57,58]};
%}


%for sb=5%:length(subs)
    %cd(['C:\shared\mong_eeg\',subs{sb}]);
    cd(['C:\shared\greg_eegfmri']);
    sets = dir('*set'); 
    for st=1:length(sets)
        EEG = pop_loadset(sets(st).name);  
        if st==1 ; merged = EEG; else merged = pop_mergeset(EEG,merged); end
    end
    
    filteeg = eegfiltfft(merged.data,EEG.srate,2,80); 
    [weights,sphere] = runica(filteeg,'maxsteps',128); 
    fullweights = weights*sphere;
    save('fullweights','fullweights'); 
    winv = pinv(weights*sphere); 
    figure,for i=1:64; subplot(5,13,i) ; topoplot(winv(:,i),EEG.chanlocs) ;title(i); end
    
    %fullweights = load('fullweights.mat') ; fullweights = fullweights.fullweights; 
    newmerged = merged; newmerged.data = fullweights*merged.data; 
    trigs = {'S 21','S 22','S 24'}; 
    %trigs = {'S  2'};
    
    clear ersp
    for i=1:length(trigs);  disp(i); 
        epica = pop_epoch(newmerged,{trigs{i}},[-2,5]); 
        for j=1:64
                [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(j,:,:)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                        'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',200) ; 
        end
    end
    figure,for i=1:64 ; subplot(5,13,i); imagesc(squeeze(mean(ersp(3,i,:,:),1)),[-5,5]); title(i); axis xy; end 
    
    %{
    [s,f] = spectopo(newmerged.data,0,newmerged.srate); 
    winv = pinv(fullweights); 
    for i=1:64 ; subplottight(10,14,i*2) ; plot(s(i,:)); text(60,0,['comp ',num2str(i)]); set(gca,'XTick',[],'YTick',[]); end 
    for i=1:64 ; subplottight(10,14,i*2-1) ; topoplot(winv(:,i),merged.chanlocs); end
    %}
    
%end
