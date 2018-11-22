clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 

% what is it that distinguishes good from bad components? correlation with
% template, what else? you can't correlate with template if its obscured by
% noise, however. 


postchans = [60,61,62,63,64,29,30,31,23,56,24,57,25,58,26,59,27,18,52,19,53,20,54,21];
badchans = zeros(1,64); badchans(postchans) = 1; badchans = find(badchans==0); 
%topoplot(zeros(1,64),eeg.chanlocs,'emarker2',{postchans,'o','b'})

cd(['E:\clean_allres\vincent']) ; ls 
allrawersp = load('allrawersp'); allrawersp = allrawersp.allrawersp; 
vmersp = squeeze(mean(mean(mean(allrawersp(:,:,postchans,:,:))))); 
timeinds = 15:90; freqinds = 5:60; 

for sb=1:length(subs); disp(sb) 
    cd(['E:\clean_allres\',subs{sb}]) ; ls 
 
    eeg = pop_loadset('interp_flip_eeg.set'); 
    cleanweights = load('cleanweights'); cleanweights = cleanweights.cleanweights; 
    weights = cleanweights{1}; sphere = cleanweights{2}; 
    winv = pinv(weights*sphere); 
    neweeg = eeg;
    acts = weights*sphere*eeg.data;     
    neweeg.data = acts; 
   
    trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'}; 
    clear ersp erspboot; 
    for tr=1:length(trigs)
    epica = pop_epoch(neweeg,{trigs{tr}},[-.85,3]); 
    [sv,si] = sort(zscore(squeeze(mean(std(epica.data,0,2),1))),'descend');
     for i=1:64
            [ersp(tr,i,:,:),itc,powbase,times,freqs,erspboot(tr,i,:,:),~] = newtimef(squeeze(epica.data(i,:,si(10:end))),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',0,'verbose','off','timesout',100,'alpha',0.05) ; 
     end
    end 
    
    
     erspmask = zeros(size(ersp)); 
     for tr=1:6
     for i=1:64
         for j=1:60
             for k=1:100
                if ersp(tr,i,j,k) < erspboot(tr,i,j,1) || ersp(tr,i,j,k) > erspboot(tr,i,j,2)
                    erspmask(tr,i,j,k) = 1; 
                end
             end
         end
     end
     end
     
     postweights = (mean(abs(winv(postchans,:)),1)./mean(abs(winv(badchans,:)),1)); 
     
     alpha = 8:14; gamma = 50:80; 
     malpha = squeeze(median(mean(mean(erspmask(:,:,6:12,times>0.5 & times<2),1),3),4)); 
     mgamma = squeeze(median(mean(mean(erspmask(:,:,25:45,times>0.5 & times<2),1),3),4)); 
     mnoise = squeeze(median(mean(mean(erspmask(:,:,25:60,times>2.1),1),3),4));
     
    mersp = squeeze(mean(ersp,1)); 
    for i=1:64 ; corrs(i) = corr2(squeeze(mersp(i,freqinds,timeinds)),vmersp(freqinds,timeinds)); end

     
    mboth = (malpha+mgamma - mnoise) .*corrs;
     
    [sv,si] = sort(mboth,'descend'); 
    figure,for i=1:64 ; subplot(8,16,i*2) ; imagesc(squeeze(mersp(si(i),:,:)),[-2,2]) ; title(si(i)); axis xy ; colormap jet; end
    for i=1:64 ; subplot(8,16,i*2-1) ; topoplot(winv(:,si(i)),eeg.chanlocs); end
    suptitle(subs{sb}); 
    
    sortcomps = si;
    save('sortcomps','sortcomps'); 
        
    
    % look for components with sustained, statistically significant
    % response to all stimulus types
    % remove bad trials and noisy components. 
    % first create the component template using raw data only (electrodes)
   
    
end




