clear all ; close all ; 

subs = {'b_alex','b_dina','b_genevieve','b_jeremie','b_karl','b_russell','b_sukhman','b_tegan','b_valerie'};

comps = {[17,22,15],[15,18,35],[17,35,5,9],[6,26,19],[10,16,15],[12,19,9],[9,13],[25,10,23],[17,7]};

for sb=1%:9
    cd(['e:/nimg_pool/',subs{sb}]); ls  
    
    eeg = pop_loadset('cleanfilt.set');
    eeg.data = eegfiltfft(eeg.data,eeg.srate,1,120); 
    %[weights,sphere] = runica(eeg.data,'maxsteps',128); 
    %winv = pinv(weights*sphere); 
    
    %filtcomps_1_100 = {}; filtcomps_1_100{1} = weights; filtcomps_1_100{2} = sphere; 
    %save('filtcomps_1_100','filtcomps_1_100'); 
    
    filtcomps_1_100 = load('filtcomps_1_100'); filtcomps_1_100 = filtcomps_1_100.filtcomps_1_100; 
    weights = filtcomps_1_100{1}; sphere = filtcomps_1_100{2}; 
    winv = pinv(weights*sphere); 
    
    stims = {'S  1','S  2','S  3'};
    neweeg = eeg; neweeg.data = weights*sphere*eeg.data; 
    clear ersp; 
    for st=1:length(stims)
    epochs = pop_epoch(neweeg,{stims{st}},[-1.5,5.5]);
    noiseepochs = pop_epoch(eeg,{stims{st}},[-2.5,6.5]);
    [sv,si] = sort(squeeze(max(max(noiseepochs.data,[],2),[],1)),'descend');    
    for i=1:64 ; disp([subs{sb},' ',num2str(i)]); 
            [ersp(st,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochs.data(i,:,si(7:end))),epochs.pnts,[epochs.xmin,epochs.xmax],epochs.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',100) ; 
    end
    
    end    
    subplot(3,3,sb) ; imagesc(squeeze(mean(mean(ersp(:,comps{sb}(1:2),:,:),1),2)),[-4,4]) ; axis xy  ;colormap jet;     
    allersp(sb,:,:,:,:) = squeeze(mean(ersp(:,comps{sb}(1:2),:,:),2)); 
       
    cd e:/nimg_pool/saved ; 
    times_9 = times; save('times_9','times_9');
    freqs_9 = freqs ; save('freqs_9','freqs_9'); 
    
    swinv = zeros(size(winv)); 
    swinv(:,comps{sb}(1:2)) = winv(:,comps{sb}(1:2)); 
    invacts = swinv*neweeg.data; 
    
    filtgamma= eegfiltfft(invacts,eeg.srate,40,80); 
    filtalpha = eegfiltfft(invacts,eeg.srate,8,25); 
    
    for st=1:length(stims)
        
    epeeg = eeg; epeeg.data = filtgamma; epgamma = pop_epoch(epeeg,{stims{st}},[-1,6]); 
    epeeg = eeg; epeeg.data = filtalpha; epalpha = pop_epoch(epeeg,{stims{st}},[-1,6]); 
    [sv,si] = sort(squeeze(mean(std(epgamma.data,0,2),1)),'descend'); 
    
    mgamma(sb,st,:) = squeeze(mean(mean(abs(epgamma.data(:,epgamma.times>0 & epgamma.times<5000,si(3:end))),2),3)) - squeeze(mean(mean(abs(epgamma.data(:,epgamma.times<0,si(3:end))),2),3));
    malpha(sb,st,:) = squeeze(mean(mean(abs(epalpha.data(:,epalpha.times>0 & epalpha.times<5000,si(3:end))),2),3)) - squeeze(mean(mean(abs(epalpha.data(:,epalpha.times<0,si(3:end))),2),3));

    end
    
    %{
    % FILTER AND RESAMPLE RAW DATA
    eeg = pop_loadbv('.',outs(1).name); 
    eeg = pop_chanedit(eeg,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
    eeg = pop_resample(eeg,250); 
    eeg.data = eeg.data - eegfiltfft(eeg.data,eeg.srate,59,61) - eegfiltfft(eeg.data,eeg.srate,84,86); 
    pop_saveset(eeg,'outside_cleanfilt.set');
    %}
    
    %{
    eeg = pop_loadset('outside_cleanfilt.set') ;  
    filtcomps = load('outside_filtcomps_1_100.mat'); weights = filtcomps.outside_filtcomps_1_100{1}; sphere = filtcomps.outside_filtcomps_1_100{2};
    winv = pinv(weights*sphere); 
    neweeg = eeg; neweeg.data = weights*sphere*eeg.data; 
    %}
    
    %{
    % CALCULATE FILTERED GAMMA AND ALPHA/BETA WITH DIFFERENT #COMPONENTS
    sortcomps_1_100 = load('outside_sortcomps_1_100'); sortcomps = sortcomps_1_100.outside_sortcomps_1_100; 
    epochs = pop_epoch(neweeg,{'S  2'},[-1,4]); 
    stds = squeeze(mean(std(epochs.data,0,2),1)); 
    [sv,si] = sort(stds,'ascend'); 
    

    swinv = zeros(size(winv)); 
    swinv(:,sortcomps(1:5)) = winv(:,sortcomps(1:5)); 
    invacts = swinv*neweeg.data; 
    neweeg.data = invacts ; 
    pop_saveset(neweeg,'filtcomps_5.set'); 
    %}
    
    % FILTER INDIVIDUAL # COMPONENTS
    %{
    for i=1:64
    
    swinv = zeros(size(winv)); 
    swinv(:,sortcomps(1:i)) = winv(:,sortcomps(1:i)); 
    invacts = swinv*neweeg.data; 
    
    filt_gamma = eegfiltfft(invacts,eeg.srate,40,60); 
    filt_alpha = eegfiltfft(invacts,eeg.srate,8,25); 
    
    eeg_gamma = eeg; eeg_gamma.data = abs(filt_gamma); 
    ep_gamma = pop_epoch(eeg_gamma,{'S  2'},[-1,4]); 
    m_gamma = squeeze(mean(ep_gamma.data(:,epochs.times>0 & epochs.times<3000,:),2)) - squeeze(mean(ep_gamma.data(:,epochs.times<0,:),2)); 
    
    eeg_alpha = eeg; eeg_alpha.data = abs(filt_alpha); 
    ep_alpha = pop_epoch(eeg_alpha,{'S  2'},[-1,4]); 
    m_alpha = squeeze(mean(ep_alpha.data(:,epochs.times>0 & epochs.times<3000,:),2)) - squeeze(mean(ep_alpha.data(:,epochs.times<0,:),2)); 

    allsub_gamma_9(sb,i,:) = squeeze(mean(m_gamma(:,si(1:end)),2)); 
    allsub_alpha_9(sb,i,:) = squeeze(mean(m_alpha(:,si(1:end)),2)); 
    disp(subs{sb}); 
    
    end
    %}
    
    
    % CORRELATE NRF WITH SINGLE TRIALS     
    %{
    eeg = pop_loadset('outside_cleanfilt.set') ;  
    filtcomps = load('outside_filtcomps_1_100.mat'); weights = filtcomps.outside_filtcomps_1_100{1}; sphere = filtcomps.outside_filtcomps_1_100{2};
    winv = pinv(weights*sphere); 
    neweeg = eeg; neweeg.data = weights*sphere*eeg.data; 
    
    epochs = pop_epoch(neweeg,{'S  2'},[-1,4]); 
    clear bersp ersp; 
    for i=1:64 ; disp([subs{sb},' ',num2str(i)]); 
        for j=1:15
            [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochs.data(i,:,j)),epochs.pnts,[epochs.xmin,epochs.xmax],epochs.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off','timesout',100) ; 
        end
    end
    bersp = ersp - repmat(mean(ersp(:,:,:,times<0),4),[1,1,1,100]);
    
    gammabersp = squeeze(mean(bersp(:,:,15:40,:),3));
    alphabersp = squeeze(mean(bersp(:,:,7:12,:),3));     
    for i=1:64
        gammacorrs(i,:) = corr(squeeze(gammabersp(i,:,:))',gamma_nrf); 
        alphacorrs(i,:) = corr(squeeze(alphabersp(i,:,:))',alpha_nrf); 
    end   
    stds = squeeze(mean(mean(std(bersp,0,4),3),1)); 
    [~,stdi] = sort(stds,'descend');    
    [sv,si] = sort(squeeze(mean((alphacorrs)+(gammacorrs),2)),'descend'); 
       
    figure,for i=1:64 ; subplottight(7,20,i*2) ; imagesc(squeeze(mean(bersp(si(i),stdi(1:end),:,:),2)),[-4,4]) ; axis xy; colormap jet; set(gca,'XTickLabel',[]','YTickLabel',[]); end 
    for i=1:64 ; subplottight(7,20,i*2-1) ; topoplot(winv(:,si(i)),eeg.chanlocs); text(-.35,.6,['component ',num2str(i)]);  end       
    outside_sortcomps_1_100 = si ; save('outside_sortcomps_1_100','outside_sortcomps_1_100'); 
    %}
    
    % RUN ICA
    %{
    filteeg = eegfiltfft(eeg.data,eeg.srate,1,100); 
    [weights,sphere] = runica(filteeg(:,1:end),'maxsteps',128); 
    winv = pinv(weights*sphere); 
 
    outside_filtcomps_1_100 = {}; outside_filtcomps_1_100{1} = weights; outside_filtcomps_1_100{2} = sphere; 
    save('outside_filtcomps_1_100','outside_filtcomps_1_100'); 
    
    neweeg = eeg; neweeg.data = weights*sphere*eeg.data; 
    epochs = pop_epoch(neweeg,{'S  2',},[-1,4]); 
    for i=1:64
            [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochs.data(i,:,:)),epochs.pnts,[epochs.xmin,epochs.xmax],epochs.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',100) ; 
    end
    figure,for i=1:64 ; subplottight(7,20,i*2) ; imagesc(squeeze(ersp(i,:,:)),[-4,4]) ; axis xy; colormap jet; set(gca,'XTickLabel',[]','YTickLabel',[]); end 
    for i=1:64 ; subplottight(7,20,i*2-1) ; topoplot(winv(:,i),eeg.chanlocs); text(-.35,.6,['component ',num2str(i)]);  end 
    %}

    
    % SAVE NRF
    %{
    eeg = pop_loadset('outside_cleanfilt.set') ;  
    filtcomps = load('outside_filtcomps_1_100.mat'); weights = filtcomps.outside_filtcomps_1_100{1}; sphere = filtcomps.outside_filtcomps_1_100{2};
    winv = pinv(weights*sphere); 
    neweeg = eeg; neweeg.data = weights*sphere*eeg.data; 
    acts = weights*sphere*eeg.data; 
    epochs = pop_epoch(neweeg,{'S  2'},[-1,4]); 
    clear ersp;
    for i=1
        for j=1:15
            [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochs.data(outside_compsm{sb}(1),:,j)),epochs.pnts,[epochs.xmin,epochs.xmax],epochs.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off','timesout',100) ; 
        end
    end
    bersp = ersp - repmat(mean(ersp(1,:,:,times<0),4),[1,1,1,100]); 
    stds = squeeze(mean(std(bersp(1,:,:,:),0,4),3)); 
    [sv,si] = sort(stds,'descend'); 
    sbcomps(sb,:,:) = squeeze(mean(bersp(1,si(2:end),:,:),2)); 
    %}
    
end
 
%cd e:/nimg_pool/saved ; 
%outside_nrf = sbcomps ; save('outside_nrf','outside_nrf'); 

%{
cd e:/nimg_pool/saved ; 
outside_allsub_gamma = allsub_gamma_9; save('outside_allsub_gamma','outside_allsub_gamma'); 
outside_allsub_alpha = allsub_alpha_9 ; save('outside_allsub_alpha','outside_allsub_alpha'); 
%}


cd e:/nimg_pool/saved ; 
mgamma_9 = mgamma; save('mgamma_9','mgamma_9');
malpha_9 = malpha; save('malpha_9','malpha_9'); 
allersp_9 = allersp; 
save('allersp_9','allersp_9'); 





