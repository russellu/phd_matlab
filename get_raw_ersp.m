clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 

% what is it that distinguishes good from bad components? correlation with
% template, what else? you can't correlate with template if its obscured by
% noise, however. 

postchans = [60,61,62,63,64,29,30,31,23,56,24,57,25,58,26,59,27,18,52,19,53,20,54,21];
ozs = [29,30,31,61,62,63,56,24,57,25,58,26,59];
badchans = zeros(1,64); badchans(postchans) = 1; badchans = find(badchans==0); 
%topoplot(zeros(1,64),eeg.chanlocs,'emarker2',{postchans,'o','b'})

for sb=1:length(subs); disp(sb) 
    cd(['E:\clean_allres\',subs{sb}]) ; ls 
 
    eeg = pop_loadset('interp_flip_eeg.set'); 
    sortcomps = load('sortcomps') ; 
    sortcomps = sortcomps.sortcomps; 
    allsortcomps(sb,:) = sortcomps; 
    eeg.data = eeg.data - eegfiltfft(eeg.data,eeg.srate,59,61); 
    cleanweights = load('cleanweights'); cleanweights = cleanweights.cleanweights; 
    weights = cleanweights{1}; sphere = cleanweights{2}; 
    winv = pinv(weights*sphere); 
    neweeg = eeg;
    acts = weights*sphere*eeg.data;     
    neweeg.data = acts; 
    %neweeg = eeg; 
    
    
    
    trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'}; 
    for tr=1:length(trigs)
    epica = pop_epoch(neweeg,{trigs{tr}},[-.85,3]); 
    [sv,si] = sort(zscore(squeeze(mean(std(epica.data,0,2),1))),'descend');
     for i=1:64
            [ersp(sb,tr,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(i,:,si(10:end))),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',0,'verbose','off','timesout',100) ; 
     end
    end 
    
end
clear subersp; 
for i=1:24 
    subplot(4,8,i) ; 
    imagesc(squeeze(mean(mean(ersp(i,[1,3,5],allsortcomps(i,1:5),:,1:end-20),2),3)),[-2,2]) ; axis xy ; colormap jet; title(subs{i}); 
    subersp(i,:,:) = squeeze(mean(mean(ersp(i,:,allsortcomps(i,1:5),:,:),2),3)); 
end



%{
allrawersp = ersp; 
save('allrawersp','allrawersp'); 
%}


