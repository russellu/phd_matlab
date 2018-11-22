clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 

badcs = {[],[],[],[3,4],[5],[],[2,3],[],[],[5],[],[],[],[],[],[3],[],[],[],[],[],[4,5],[],[]};

for sb=1%:length(subs)
    disp(sb); 
    cd(['E:\clean_allres\',subs{sb}]) ; ls 

    eeg = pop_loadset('interp_flip_eeg.set'); 
    eeg.data = eeg.data - eegfiltfft(eeg.data,eeg.srate,59,61); 
    cleanweights = load('cleanweights'); cleanweights = cleanweights.cleanweights; 
    weights = cleanweights{1}; sphere = cleanweights{2}; 
    winv = pinv(weights*sphere); 
    neweeg = eeg;
    acts = weights*sphere*eeg.data;     
    sortcomps = load('sortcomps') ; 
    sortcomps = sortcomps.sortcomps; 
    eeg.data = acts; 
    
    trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'}; 
    for tr=1:length(trigs)
    epica = pop_epoch(eeg,{trigs{tr}},[-.5,2.5]); 
    [sv,si] = sort(zscore(squeeze(mean(std(epica.data,0,2),1))),'descend');
    epica.data(:,:,si(1:10)) = []; 
    for i=1:5
        [ersp(sb,tr,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(sortcomps(i),:,:)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',0,'verbose','off','timesout',100) ; 
    end
    end
end

for i=1:5 ; figure; for j=1:24; subplot(5,5,j) ; imagesc(squeeze(mean(ersp(j,:,i,:,:),2)),[-1,1]) ; axis xy; colormap jet; title(j); end; end

clear mersp; 
for i=1:24
   badsi = badcs{i}; 
   goods = zeros(1,5); 
   goods(badsi) = 1; goods = find(goods==0); 
   mersp(i,:,:,:) = squeeze(mean(ersp(i,:,goods,:,:),3)); 
end

submersp = mersp; save('submersp','submersp'); 


