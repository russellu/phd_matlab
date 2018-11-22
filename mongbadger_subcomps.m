clear all ; close all; 

subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'} ;
comps = {[4,5,6,7,8,9,10,11,12,13,14,15,17,20,21,22,23,25,26,28,30,38],[5,7,8,9,13,17,18,19,20,25,26,28,29,37,47,58],...
         [4,5,6,7,8,9,11,13,14,15,16,17,20,21,23,25,26,27,28,35],[3,5,8,9,10,12,17,18,19,21,22,27,29,36],...
         [4,8,9,11,12,13,14,16,20,21,25,27,36],[5,6,7,8,9,10,12,13,14,16,20,26,35,37],...
         [1,6,8,9,11,13,15,17,24,25,28,33,35,39],[3,5,7,8,10,15,16,20,21,23,26,27,28,29,30,38,40],...
         [3,6,8,10,11,16,19,20,22,25,30,36,42]};
%{
subs = {'MONG_01_RB','MONG_02_DP','MONG_03_CG','MONG_05_SG','MONG_06_TS'}; 
comps = {[3,5,6,7,8,10,11,14,16,17,21,22,24,25,26,30,31],[2,4,9,10,11,15,16,17,18,20,22,27,29,42,48],...
         [2,4,7,9,11,12,15,16,18,23,24,29,30,32,33,34,36,39,40,41,46],[5,6,7,8,9,10,11,15,16,17,18,19,20,24,25,37,48,49],...
         [6,10,13,15,16,18,19,21,23,25,28,41,43,52,57,58]};
%}


for sb=1:length(subs)
    cd(['C:\shared\badger_eeg\',subs{sb}]);
    sets = dir('*gamma*set'); 
    for st=1:length(sets)
        EEG = pop_loadset(sets(st).name);  
        if st==1 ; merged = EEG; else merged = pop_mergeset(EEG,merged); end
    end
    %{
    filteeg = eegfiltfft(merged.data,EEG.srate,3,128); 
    [weights,sphere] = runica(filteeg,'maxsteps',128); 
    fullweights = weights*sphere;
    save('fullweights','fullweights'); 
    winv = pinv(weights*sphere); 
    figure,for i=1:64; subplottight(5,13,i) ; topoplot(winv(:,i),EEG.chanlocs) ; end
    %}
    fullweights = load('fullweights.mat') ; fullweights = fullweights.fullweights; 
    %newmerged = merged; newmerged.data = fullweights*merged.data; 
    acts = fullweights*merged.data; goods = comps{sb}; 
    bads = zeros(1,64); bads(goods) = 1 ; bads = find(bads==0); 
    acts(bads,:) = 0; 
    invdat = pinv(fullweights)*acts; 
    newmerged = merged; newmerged.data = invdat; 
    trigs = {'S  1','S  2','S  3'}; 
    %trigs = {'S  2'};  
    clear ersp
    for i=1:length(trigs);  disp(i); 
        epica = pop_epoch(newmerged,{trigs{i}},[-1.5,6]); 
        for j=1:64
                [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(j,:,:)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                        'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',200) ; 
        end
    end
    %figure,for i=1:64 ; subplot(6,13,i); imagesc(squeeze(mean(ersp(1,i,:,:),1)),[-5,5]); axis xy; end 
    figure,topoplot(squeeze(mean(mean(mean(ersp(1:2,:,5:12,times>0 & times<5),3),4),1)),merged.chanlocs);
end
