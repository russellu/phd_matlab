clear all ; close all ; 
subs = {'lowcont_carmen','lowcont_lyndis','lowcont_mael','lowcont_lorne','lowcont_lyes','lowcont_greg','lowcont_angelina','lowcont_cesar','lowcont_charles','lowcont_alex','lowcont_janie','lowcont_russell','lowcont_anne'} ; 
stims = {'S  2','S  4','S 12','S 14','S 22','S 24','S 32','S 34','S 42','S 44'} ;
bades = {[7,32,15,39],[4,33,5,38],[12,4,46,37,16,19,59,62,6],[15,10,20],[10],[19,32,62,18,56],[32,1,15],[62,12,1,40,19,4],[46,37,12,58,19],[19,58,32,16,37,18],[32,19,62,58,46],[62,17,59,28],[6,62,31,19,30,64]}; 
firstpass = {[5,6,18],[10,11],[8,14,24],[5],[13,22],[6,11],[11,22,49],[18,32,45],[12,21],[11,23],[8],[7,15,22],[9,17]};

cd e:/saved ;
mean_nrf = load('mean_nrf'); mean_nrf = mean_nrf.mean_nrf; 
gamma_ts = squeeze(mean(mean(mean_nrf(:,:,18:35,:),1),3)); 
alpha_ts = squeeze(mean(mean(mean_nrf(:,:,5:8,:),1),3)); 

for sb=1:length(subs)
    
    cd(['E:/jly_orientation/',subs{sb}]) ;  
    %{
    ls *vhdr    
    vhdrs = dir('*vhdr');    
    for vhdr=1:length(vhdrs)
       if vhdr==1; merged = pop_loadbv('.',vhdrs(vhdr).name); else merged = pop_mergeset(merged,pop_loadbv('.',vhdrs(vhdr).name)); end
    end
    merged = pop_chanedit(merged,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;    
    % isolate bad channels
    resmerged = pop_resample(merged,250); 
    resmerged.data = resmerged.data - eegfiltfft(resmerged.data,resmerged.srate,59,61) - eegfiltfft(resmerged.data,resmerged.srate,84,86) ;  
    %}
    %{
    diffs = abs(diff(merged.data,1,2)); 
    skews = zscore(skewness(diffs,0,2)); 
    sumdiffs = zscore(sum(diffs,2)); 
    [sv,si] = sort(sumdiffs + skews,'descend'); 
    figure,subplot(2,1,1) ; bar(sv) ; set(gca,'XTick',1:64,'XTickLabel',si); 
    subplot(2,1,2) ; topoplot(sumdiffs+skews,merged.chanlocs,'electrodes','numbers'); suptitle(subs{sb}); 
    %}
    
    %{
    resmerged = pop_interp(resmerged,bades{sb},'spherical');
    pop_saveset(resmerged,'merged2.set');         
    eeg=pop_loadset('merged2.set'); 
    mergefilt = eeg ; mergefilt.data = eegfiltfft(eeg.data,eeg.srate,1,80)  ; 
    [weights,sphere] = runica(mergefilt.data(:,1:end),'maxsteps',256,'stop',1e-9);     
    eegcomps2{1} = weights; eegcomps2{2} = sphere; save('eegcomps2','eegcomps2'); 
    %}   
    
    eeg=pop_loadset('merged2.set'); 
    eegweights = load('eegcomps2'); eegweights = eegweights.eegcomps2; 
    newmerged = eeg ; newmerged.data = eegweights{1}*eegweights{2}*eeg.data ; 
    clear allersp; 
    for s=1:length(stims)
        allep = pop_epoch(newmerged,{stims{s}},[-2,20]) ; disp(s); 
        for i=1:64 
            for j=1:size(allep.data,3)
                [allersp(s,i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,j)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                        'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',NaN,'verbose','off','timesout',200) ; 
            end    

        end
    end
    bersp = allersp - repmat(mean(allersp(:,:,:,:,times<0),5),[1,1,1,1,200]);
    
    mgammabersp = squeeze(mean(bersp(:,:,:,18:35,:),4));
    malphabersp = squeeze(mean(bersp(:,:,:,5:8,:),4)); 
    
    clear gammacorrs alphacorrs
    for i=1:size(mgammabersp,1)
        for j=1:size(mgammabersp,2)
            for k=1:size(mgammabersp,3)
                gammacorrs(i,j,k) = corr2(squeeze(mgammabersp(i,j,k,:)),gamma_ts(i,:)'); 
                alphacorrs(i,j,k) = corr2(squeeze(malphabersp(i,j,k,:)),alpha_ts(i,:)'); 
            end
        end
    end
    gammacorrs(isnan(gammacorrs)) = 0;
    alphacorrs(isnan(alphacorrs)) = 0; 
    mcorrs = squeeze(mean(mean(gammacorrs,1),3)) + squeeze(mean(mean(alphacorrs,1),3));    
    [sv,si] = sort(mcorrs,'descend'); 
    mbersp = squeeze(mean(mean(bersp,1),3));     
    newcomps_si = si ; save('newcomps_si','newcomps_si');     
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mbersp(si(i),:,:)),[-3,3]) ; colormap jet; axis xy;  end
        
    %{
    clear allersp ; 
    newmerged = eeg ; newmerged.data = eegweights{1}*eegweights{2}*eeg.data ;  clear allersp; 
    for s=1:length(stims)
    allep = pop_epoch(newmerged,{stims{s}},[-2,20]) ;
    for i=1:64 
            [allersp(s,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',0,'verbose','off','timesout',200) ; 
    end    
    end
    %figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(allersp(:,i,:,:)),1),[-5,5]) ; colormap jet ;  title(i) ; axis xy ; end ; suptitle(subs{sb}); 
    mvals(sb,:,:,:,:) = allersp; 
    %}
end
%{
for i=1:13 ; figure; for j=1:64 ;subplot(5,13,j) ; imagesc(squeeze(mean(mvals(i,:,j,:,:),2)),[-3,3]) ; axis xy ; colormap jet; title(j); end; suptitle(subs{i}); end
for i=1:13 ; mean_nrf(i,:,:,:) = squeeze(mean(mvals(i,:,firstpass{i},:,:),3)) ; end
cd e:/saved ; save('mean_nrf','mean_nrf'); 
%}






