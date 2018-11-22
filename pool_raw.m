clear all ; close all ; 
subs = {'a_alex','a_charest','a_esteban','a_fabio','a_gab','a_gabriella','a_genevieve','a_gina','a_guillaume','a_jeremie','a_julie','a_katrine','a_lisa','a_marc',...
    'a_marie','a_mathieu','a_maxime','a_mingham','a_patricia','a_po','a_russell','a_sunachakan','a_tah','a_vincent'};

postelecs = [29,30,31,60,61,62,63,64,23,56,24,57,25,58,26,59,27];

compsm = {[4,5,7,8],[6,2,5],[7,6,3],[6,9,14],[5,7,11,12],[5,8,10],[9,6,3],[3,9,6],[1,6,19],[3,15,4],[6,10,12,11],[4,2,10,18],[3,12,11],[5,7,12],[5,4,12],[1,15,11],[2,4,8],[10,3,7,5],[2,5,8],[4,9,8],[20,15,3,21],[5,1,24],[8,5,6],[2,4,7]};
cd E:\nimg_pool\saved ; nrf = load('nrf'); nrf = nrf.nrf; 
gamma_nrf = squeeze(mean(mean(nrf(:,15:40,:),1),2)); alpha_nrf = squeeze(mean(mean(nrf(:,7:12,:),1),2)); 
for sb=1:24
    cd(['e:/nimg_pool/',subs{sb}]); ls  
    
    eeg = pop_loadset('cleanfilt.set') ; 
    
    eeg.data =eegfiltfft(eeg.data,eeg.srate,7,100); 
    
    
    filtcomps = load('filtcomps_1_100.mat'); weights = filtcomps.filtcomps_1_100{1}; sphere = filtcomps.filtcomps_1_100{2};
    winv = pinv(weights*sphere); 
    neweeg = eeg; neweeg.data = weights*sphere*eeg.data; 
    
    sortcomps_1_100 = load('sortcomps_1_100'); sortcomps = sortcomps_1_100.sortcomps_1_100; 
    epochs = pop_epoch(neweeg,{'S 11','S 12','S 13','S 14','S 15','S 16'},[-.85,2.85]); 
    stds = squeeze(mean(std(epochs.data,0,2),1)); 
    [sv,si] = sort(stds,'ascend'); 
    
    
    
    % GET VARIANCE ACCOUNTED FOR BY EACH COMPONENT    
    %clear vars pvaf; 
    for i=1:64 ; disp(i); 
        swinv = zeros(size(winv)); 
    
        swinv(:,sortcomps(1:i)) = winv(:,sortcomps(1:i)); 
        invacts = swinv*neweeg.data; 
        
        %vars(i) = mean(var(eeg.data-invacts))/mean(var(eeg.data));
        
        %var1 = sum(sum(abs(diff(eeg.data - invacts,1,2)),1),2);
        %var2 = sum(sum(abs(diff(eeg.data,1,2)),1),2);
        
        %rat = 100*(var1/var2);
        %pvaf(i) = rat; 
        pvaf(sb,i) = 100-100*mean(var(eeg.data-invacts))/mean(var(eeg.data));
      
    end
    
    %pvaf(comp) = 100-100*mean(var(data - back_proj))/mean(var(data));
    
    
    
    % SAVE THE TOP 5 COMPONENTS AS RAW DATASET
    %{
    swinv = zeros(size(winv)); 
    swinv(:,sortcomps(1:5)) = winv(:,sortcomps(1:5)); 
    invacts = swinv*neweeg.data; 
    neweeg.data = invacts;         
    pop_saveset(neweeg,'filtcomps_5.set');  
    %}
    
    
    
    
    %{
    trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'};
    clear ersp;
    for tr=1:length(trigs)
        sortcomps_1_100 = load('sortcomps_1_100'); sortcomps = sortcomps_1_100.sortcomps_1_100; 
        epochs = pop_epoch(neweeg,{trigs{tr}},[-.85,2.85]); 
        stds = squeeze(mean(std(epochs.data,0,2),1)); 
        [sv,si] = sort(stds,'ascend'); 
        for i=1:64 ; disp([subs{sb},' ',num2str(i)]); 
                [ersp(i,tr,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochs.data(sortcomps(i),:,si(1:end-10))),epochs.pnts,[epochs.xmin,epochs.xmax],epochs.srate,0,...
                        'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',100) ; 
        end   
    end
    allsub_stims(sb,:,:,:,:) = ersp; 
    %}
    
    %cd e:/nimg_pool/saved ; save('times','times') ; save('freqs','freqs'); 
    
    %{
    % visualize sorted and unsorted components 
    
    sortcomps_1_100 = load('sortcomps_1_100'); sortcomps = sortcomps_1_100.sortcomps_1_100; 
    epochs = pop_epoch(neweeg,{'S 11','S 12','S 13','S 14','S 15','S 16'},[-.85,2.85]); 
    stds = squeeze(mean(std(epochs.data,0,2),1)); 
    [sv,si] = sort(stds,'ascend'); 
    
    clear ersp;
    for i=1:64 ; disp([subs{sb},' ',num2str(i)]); 
            [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochs.data(sortcomps(i),:,si(1:end-100))),epochs.pnts,[epochs.xmin,epochs.xmax],epochs.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',100) ; 
    end   
    figure,for i=1:64 ; subplottight(4,32,i*2) ; imagesc(squeeze(ersp(i,:,:)),[-4,4]) ; axis xy; colormap parula; set(gca,'XTickLabel',[]','YTickLabel',[]); end 
    for i=1:64 ; subplottight(4,32,i*2-1) ; topoplot(winv(:,sortcomps(i)),eeg.chanlocs); text(0.2,.55,['component ',num2str(i)]); colormap parula;  end 
    %}
    
    %{
    % STIMULUS SPECIFIC FILTERED POWER
    sortcomps_1_100 = load('sortcomps_1_100'); sortcomps = sortcomps_1_100.sortcomps_1_100; 
    trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'}; 
    for tr=1:length(trigs)
    
    epochs = pop_epoch(neweeg,{trigs{tr}},[-.85,2.85]); 
    stds = squeeze(mean(std(epochs.data,0,2),1)); 
    [sv,si] = sort(stds,'ascend'); 
    
    
    swinv = zeros(size(winv)); 
    swinv(:,sortcomps(1:5)) = winv(:,sortcomps(1:5)); 
    invacts = swinv*neweeg.data; 
    
    filt_gamma = eegfiltfft(invacts,eeg.srate,40,100); 
    filt_alpha = eegfiltfft(invacts,eeg.srate,8,25); 
    
    eeg_gamma = eeg; eeg_gamma.data = abs(filt_gamma); 
    ep_gamma = pop_epoch(eeg_gamma,{trigs{tr}},[-.85,2.85]); 
    m_gamma = squeeze(mean(ep_gamma.data(:,epochs.times>0 & epochs.times<2000,:),2)) - squeeze(mean(ep_gamma.data(:,epochs.times<0,:),2)); 
    
    eeg_alpha = eeg; eeg_alpha.data = abs(filt_alpha); 
    ep_alpha = pop_epoch(eeg_alpha,{trigs{tr}},[-.85,2.85]); 
    m_alpha = squeeze(mean(ep_alpha.data(:,epochs.times>0 & epochs.times<2000,:),2)) - squeeze(mean(ep_alpha.data(:,epochs.times<0,:),2)); 

    allsub_allstim_gamma(sb,tr,:) = squeeze(mean(m_gamma(:,si(1:end-15)),2)); 
    allsub_allstim_alpha(sb,tr,:) = squeeze(mean(m_alpha(:,si(1:end-15)),2)); 
    disp(subs{sb});  
    end
    %}
   
    %{
    % CALCULATE FILTERED GAMMA AND ALPHA/BETA WITH DIFFERENT #COMPONENTS
    sortcomps_1_100 = load('sortcomps_1_100'); sortcomps = sortcomps_1_100.sortcomps_1_100; 
    epochs = pop_epoch(neweeg,{'S 11','S 12','S 13','S 14','S 15','S 16'},[-.85,2.85]); 
    stds = squeeze(mean(std(epochs.data,0,2),1)); 
    [sv,si] = sort(stds,'ascend'); 
    
    for i=1:64
    
    swinv = zeros(size(winv)); 
    swinv(:,sortcomps(1:i)) = winv(:,sortcomps(1:i)); 
    invacts = swinv*neweeg.data; 
    
    filt_gamma = eegfiltfft(invacts,eeg.srate,40,100); 
    filt_alpha = eegfiltfft(invacts,eeg.srate,8,25); 
    
    eeg_gamma = eeg; eeg_gamma.data = abs(filt_gamma); 
    ep_gamma = pop_epoch(eeg_gamma,{'S 11','S 12','S 13','S 14','S 15','S 16'},[-.85,2.85]); 
    m_gamma = squeeze(mean(ep_gamma.data(:,epochs.times>0 & epochs.times<2000,:),2)) - squeeze(mean(ep_gamma.data(:,epochs.times<0,:),2)); 
    
    eeg_alpha = eeg; eeg_alpha.data = abs(filt_alpha); 
    ep_alpha = pop_epoch(eeg_alpha,{'S 11','S 12','S 13','S 14','S 15','S 16'},[-.85,2.85]); 
    m_alpha = squeeze(mean(ep_alpha.data(:,epochs.times>0 & epochs.times<2000,:),2)) - squeeze(mean(ep_alpha.data(:,epochs.times<0,:),2)); 

    allsub_gamma(sb,i,:) = squeeze(mean(m_gamma(:,si(1:end-100)),2)); 
    allsub_alpha(sb,i,:) = squeeze(mean(m_alpha(:,si(1:end-100)),2)); 
    disp(subs{sb}); 
    
    end
    %}
    
    %{
    % GET ALL SORTED COMPONENT ERSP 
    clear ersp;
    for i=1:64 ; disp([subs{sb},' ',num2str(i)]); 
            [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochs.data(sortcomps(i),:,si(1:end-100))),epochs.pnts,[epochs.xmin,epochs.xmax],epochs.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',100) ; 
    end   
    sbgcomps(sb,:,:,:) = ersp; 
    %}
    
    %{
    % CORRELATE TIME SERIES
    ideal = zeros(1,size(neweeg.data,2)); 
    boundarypts = 0;
    for i=1:length(types)
        if strcmpi(types{i},'boundary')
            boundarypts = boundarypts + lats(i);            
        end
        for j=1:length(trigs)
            if strcmpi(types{i},trigs{j})
                ideal(boundarypts+lats(i):boundarypts+lats(i)+eeg.srate*2) = 1; 
            end
        end
    end
    %}
    
    % CORRELATE NRF WITH SINGLE TRIALS      
    %{
    epochs = pop_epoch(neweeg,{'S 11','S 13','S 14','S 15'},[-.85,2.85]); 
    clear bersp ersp; 
    for i=1:64 ; disp([subs{sb},' ',num2str(i)]); 
        for j=1:405
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
       
    figure,for i=1:64 ; subplottight(7,20,i*2) ; imagesc(squeeze(mean(bersp(si(i),stdi(1:end-50),:,:),2)),[-4,4]) ; axis xy; colormap jet; set(gca,'XTickLabel',[]','YTickLabel',[]); end 
    for i=1:64 ; subplottight(7,20,i*2-1) ; topoplot(winv(:,si(i)),eeg.chanlocs); text(-.35,.6,['component ',num2str(i)]);  end       
    sortcomps_1_100 = si ; save('sortcomps_1_100','sortcomps_1_100'); 
    %}
    
    %{
    % RUN ICA
    filteeg = eegfiltfft(eeg.data,eeg.srate,50,80); 
    [weights,sphere] = runica(filteeg(:,1:3:end),'maxsteps',128); 
    winv = pinv(weights*sphere); 
    acts = weights*sphere*eeg.data;   
    %}   
    %{
    % SAVE NRF
    filtcomps = load('filtcomps_1_100.mat'); weights = filtcomps.filtcomps_1_100{1}; sphere = filtcomps.filtcomps_1_100{2};
    winv = pinv(weights*sphere); 
    neweeg = eeg; neweeg.data = weights*sphere*eeg.data; 
    acts = weights*sphere*eeg.data; 
    epochs = pop_epoch(neweeg,{'S 11'},[-.85,2.85]); 
    for i=1
            [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochs.data(compsm{sb}(1),:,:)),epochs.pnts,[epochs.xmin,epochs.xmax],epochs.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',100) ; 
    end
    sbcomps(sb,:,:) = ersp; 
    %}
    %{
    % VISUALIZE ALL COMPONENTS
    filtcomps = load('filtcomps_1_100.mat'); weights = filtcomps.filtcomps_1_100{1}; sphere = filtcomps.filtcomps_1_100{2};
    winv = pinv(weights*sphere); 
    neweeg = eeg; neweeg.data = weights*sphere*eeg.data; 
    acts = weights*sphere*eeg.data; 
    epochs = pop_epoch(neweeg,{'S 11'},[-.85,2.85]); 
    for i=1:64
            [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochs.data(i,:,:)),epochs.pnts,[epochs.xmin,epochs.xmax],epochs.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',100) ; 
    end
    figure,for i=1:64 ; subplottight(7,20,i*2) ; imagesc(squeeze(ersp(i,:,:)),[-4,4]) ; axis xy; colormap jet; set(gca,'XTickLabel',[]','YTickLabel',[]); end 
    for i=1:64 ; subplottight(7,20,i*2-1) ; topoplot(winv(:,i),eeg.chanlocs); text(-.35,.6,['component ',num2str(i)]);  end 
    %}
    
    %{
    % ABSOLUTE GAMMA RAW
    filteeg = eegfiltfft(eeg.data,eeg.srate,40,80); 
    filteeg = abs(filteeg); 
    epochs = pop_epoch(eeg,{'S 11'},[-.85,2.5]); 

    stdepochs = squeeze(mean(std(epochs.data(:,:,:),0,2),1)); 
    [sv,si] = sort(stdepochs,'ascend'); 
    
    m_epochs = squeeze(mean(epochs.data(:,:,si(1:end-15)),3)); 
    gamma_diff = mean(m_epochs(:,epochs.times>0 & epochs.times<2000),2) - mean(m_epochs(:,epochs.times<0),2);
    
    all_gamma_diff(sb) = mean(gamma_diff(postelecs)); 
    %}
    
    %{    
    % ERSP RAW
    %subplot(3,8,sb) ; topoplot(mean(abs(diff(eeg.data(:,:),1,2)),2),eeg.chanlocs); 
    epochs = pop_epoch(eeg,{'S 11'},[-.85,2.5]); 
    for i=1:length(postelecs)
        for j=1:135
            [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochs.data(postelecs(i),:,j)),epochs.pnts,[epochs.xmin,epochs.xmax],epochs.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off','timesout',100) ; 
        end
    end
    bersp = ersp - repmat(mean(ersp(:,:,:,times<0),4),[1,1,1,100]);
    zs = (zscore(squeeze(mean(mean(std(ersp(:,:,:,:),0,4),3),1)))); 
    [sv,si] = sort(zs,'ascend'); 

    allersp(sb,:,:,:) = squeeze(mean(bersp(:,si(1:end-15),:,:),2));
    %}
    
end
%{
% SAVE ERSP RAW
cd E:\nimg_pool\saved
a_allersp = allersp ; save('a_allersp','a_allersp'); 
%}
%{
% SAVE ABSOLUTE GAMMA RAW
cd E:\nimg_pool\saved
raw_gamma_diff = all_gamma_diff;
save('raw_gamma_diff','raw_gamma_diff'); 
%}
%{
% SAVE NRF
for i=1:24 ; subplot(4,6,i) ; imagesc(times,freqs,squeeze(sbcomps(i,:,:)),[-4,4]) ; axis xy ; colormap jet ; title(['subject ',num2str(i)]); if i==1; xlabel('time(s)'); ylabel('frequency(hz)'); end;  end
subplot(1,2,1); 
imagesc(times,freqs,squeeze(mean(sbcomps,1)),[-3,3]) ; axis xy ; colormap jet; xlabel('time(s)') ; ylabel('frequency(hz)'); 
subplot(1,2,2); 
plot(times,squeeze(mean(mean(sbcomps(:,20:40,:),1),2)),'r','LineWidth',3);  hold on ; plot(times,squeeze(mean(mean(sbcomps(:,5:10,:),1),2)),'b','LineWidth',3) ; hline(0,'k'); vline([0,2],'k'); 
legend({'gamma','alpha/beta'}); xlabel('time(s)'); ylabel('db'); xlim([-.500,2.500]); 
cd E:\nimg_pool\saved
nrf = sbcomps ; save('nrf','nrf'); 
%}

% SAVE SORTED ERSP
%{
cd E:\nimg_pool\saved
top10_n24 = sbgcomps ; save('top10_n24','top10_n24'); 
%}

% SAVE BACK PROJECTED DATA WITH DIFFERENT #COMPONENTS
%cd E:\nimg_pool\saved
%save('allsub_allstim_gamma','allsub_allstim_gamma'); 
%save('allsub_allstim_alpha','allsub_allstim_alpha'); 

% save stimulus specific ERSP
%allstim_ersp = allsub_stims; save('allstim_ersp','allstim_ersp'); 
fig = figure;
subplot(4,8,1);
imagesc(times,freqs,squeeze(mean(nrf,1)),[-2,2]) ; axis xy ; colormap parula; xlabel('time(s)'); ylabel('frequency(hz)'); 
subplot(4,8,2); 
shadedErrorBar(times,squeeze(mean(mean(nrf(:,20:end,:),1),2)),squeeze(std(mean(nrf(:,20:end,:),2),0,1))/sqrt(24),{'r'});
hold on ; 
shadedErrorBar(times,squeeze(mean(mean(nrf(:,4:12,:),1),2)),squeeze(std(mean(nrf(:,4:12,:),2),0,1))/sqrt(24),{'b'}); hline(0,'k'); 
ylim([-2.6,1.25]); xlabel('time(s)'); ylabel('power (db)'); xlim([-.7,2.7]); 
subplot(4,8,3); plot(1,'r') ; hold on ; plot(2,'b'); legend({'gamma (40-100Hz)','alpha/beta (8-25Hz)'});
set(fig,'Units','normalized');
set(fig,'Position',[0 0 1.2 .8]);

fig=figure;
subplot(4,8,1); 
bar(pvaf(:,5)); title(['mean=',num2str(round(mean(pvaf(:,5)))),'% std=',num2str(round(std(pvaf(:,5)))),'%']);
xlabel('subject'); ylabel('top 5 components %var'); 
set(fig,'Units','normalized');
set(fig,'Position',[0 0 1.2 .8]);







