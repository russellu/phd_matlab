clear all ; close all 
cd c:/shared/allres/lisa ; 
resamps = dir('resamp_vis*set') ; 
interpChans = 0 ; 
for res=1:length(resamps) ; 
    % merge the sets
    
    EEG = pop_loadset(resamps(res).name) ;
    figure,
    subplot(2,2,1) ; topoplot(double(zscore(sum(diff(EEG.data,1,2).^2,2))),EEG.chanlocs) ; 
    subplot(2,2,2) ; bar(zscore((sum(diff(EEG.data,1,2).^2,2)))) ;
    if interpChans
        chandiffs = zscore(sum(diff(EEG.data,1,2).^2,2)) ; 
        zchans = find(chandiffs>.5) ;
        EEG = pop_interp(EEG,zchans,'spherical') ; 
    end
    % concatenates
    if res ~= 1 
        merged = pop_mergeset(EEG,merged) ;
    else merged = EEG ;
    end
    

    %}
    
end
    
    resmerged = pop_resample(merged,300) ; 
    copy = resmerged ; 
    resmerged.data = eegfiltfft(resmerged.data,resmerged.srate,1,150) ; 
    resmerged = pop_runica(resmerged,'runica') ;     
    
    % apply the weights
    copy.icaact = icaact(copy.data,resmerged.icawinv,mean(copy.data,1));
    copy.icawinv = resmerged.icawinv ; copy.icasphere = resmerged.icasphere ; copy.icaweights = resmerged.icaweights ; copy.icachansind = resmerged.icachansind ; 
    eps = pop_epoch(copy,{'S 11'},[-.85,2.85]) ; 
    clear ersp ; 
    for i=1:size(eps.data,1)
        [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(eps.icaact(i,:,:)),size(eps.data,2),[eps.xmin,eps.xmax],eps.srate,0,'plotersp','off','plotitc','off',...
                                                    'baseline',0,'freqs',[1,150],'nfreqs',75,'timesout',200,'winsize',round(eps.srate/3)) ;     
    end
    
    for i=1:size(ersp,1)
        subplot(5,13,i) ; 
        imagesc(squeeze(ersp(i,:,:)),[-3,3]) ; title(num2str(i)) ; 
    end
    
   a =  [21,23,55,64] ; 
  imagesc(squeeze(mean(ersp(a,:,:),1)),[-2,2]) ; 
   
   