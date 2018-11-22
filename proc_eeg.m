cd('c:/shared/retno') ; ls 
sounds=dir('r_*vhdr') ;
for i=1:max(size(sounds)) ; 
   EEG = pop_loadbv('.',sounds(i).name) ; 
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
   EEG = pop_resample(EEG,256) ; 
   if i==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end
    
end
zsums = zscore(sum(diff(merged.data,2).^2,2)) ; 
merged = pop_interp(merged,find(zsums>1),'spherical') ; 
ica = merged ; ica.data = eegfiltfft(ica.data,ica.srate,1,128) ;
ica = pop_runica(ica,'runica') ; 
trigs = {'S 10','S 11','S 12','S 13','S 14','S 15'} ;
clear ersp ;
for i=1:size(trigs,2) 
    ep = pop_epoch(ica,{trigs{i}},[-.85,2.85]) ;
    for j=1:64 ; 
        [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(j,:,:)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,'plotersp','off','plotitc','off',...
            'freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',0) ; 
        
    end
end
figure,for i=1:64 ; subplot(5,13,i) ; topoplot(squeeze(ica.icawinv(:,i)),ica.chanlocs) ; title(i) ; end
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(ersp(:,i,:,:),1)),[-4,4]) ; title(i) ; end
for i=1:64 ; subplot(8,8,i) ; bar(zscore(ica.icawinv(:,i))) ; hline(6) ; if max(zscore(ica.icawinv(:,i))) > 7.5 ; bads(i) = 1 ; else bads(i) = 0 ; end ; end 
%goodcomps = [4,11,20,49] ; bads = ones(1,64) ; bads(goodcomps) = 0 ; 

cleaned = pop_subcomp(ica,find(bads==1),0) ; 
%{
% process the erp electrodes
for i=1:length(trigs)
    epcleans{i} = pop_epoch(cleaned,{trigs{i}},[-.85,2.85]) ; 
    epnoises{i} = pop_epoch(ica,{trigs{i}},[-.85,2.85]) ; 
    
end

t = epcleans{2} ; 
t2 = epnoises{2} ; 
tstart = find(abs(t.times)==min(abs(t.times))) ; tend = find(abs(t.times-2000)==min(abs(t.times-2000))) ; 
plot(squeeze(t.data(62,:,1))) ; hold on ; plot(smooth(squeeze(t.data(62,:,1)),10),'r') ; vline([tstart,tend]) ; 
set(gca,'XTick',1:10:length(t.times),'XTickLabel',round(t.times(1:10:end))) 
peak1 = find(t.times>80 & t.times<180) ; peak2 = find(t.times>160 & t.times<220) ; peak3 = find(t.times> 210 & t.times < 280) ; % time intervals for the 3 erp peaks (max,min,max)
figure,
for i=1:25 %size(t.data,3) ;
    erpi = squeeze(t.data(62,:,i)) ; 
    smtherp = smooth(erpi,10) ; 
    pk1 = find(smtherp(peak1)==max(smtherp(peak1))) + min(peak1) ; 
    pk2 = find(smtherp(peak2)==min(smtherp(peak2))) + min(peak2) ; 
    pk3 = find(smtherp(peak3)==max(smtherp(peak3))) + min(peak3) ; 
    allpeaks(i,:) = [pk1,pk2,pk3] ; 
    subplot(5,5,i),topoplot(double(squeeze(mean(t2.data(:,pk1-5:pk1+5,i),2))),t.chanlocs,'maplimits',[-50,50]) ; 
    % figure ; plot(erpi) ; vline([pk1,pk2,pk3]) ; 
end
postelecs = [60:64,29:31,23:27,56:59] ; 
% plot the mean erp (cleaned and non-cleaned)
for i=1:length(epcleans)
    cleani = epcleans{i}.data ; 
    noisei = epnoises{i}.data ; 
    mclean = squeeze(mean(cleani,3)) ; mnoise = squeeze(mean(noisei,3)) ; 
    pk1 = find(mclean(peak1)==max(mclean(peak1))) + min(peak1) ; 
    pk2 = find(mclean(peak2)==min(mclean(peak2))) + min(peak2) ; 
    pk3 = find(mclean(peak3)==max(mclean(peak3))) + min(peak3) ; 
    subplot(2,6,i) ; topoplot(double(squeeze(mean(cleani(:,pk1-5:pk1+5,i),2)))-double(squeeze(mean(cleani(:,pk2-5:pk2+5,i),2))),t.chanlocs,'maplimits',[-20,20],'electrodes','off') ; title(snames{i}) ; 
    subplot(2,6,i+6) ; topoplot(double(squeeze(mean(noisei(:,pk1-5:pk1+5,i),2)))-double(squeeze(mean(cleani(:,pk2-5:pk2+5,i),2))),t.chanlocs,'maplimits',[-20,20],'electrodes','off') ; title(snames{i}) ;
    allcleanerps(i,:) = squeeze(mean(mean(cleani(postelecs,:,:),1),3)) ; 
    allnoiseerps(i,:) = squeeze(mean(mean(noisei(postelecs,:,:),1),3)) ; 
end ; suptitle('denoised (top) and non-denoised topoplots peak1-peak2') ; 
plot(allcleanerps','LineWidth',3) ; legend(snames) ; 
plot(allnoiseerps','LineWidth',3) ; legend(snames) ; 
%}
%%% get the bad trial indices:
for i=1:length(epnoises)
    trialsi = squeeze(mean(epnoises{i}.data(postelecs,:,:),1)) ; 
    figure,bar(zscore(std(trialsi,0,1))) ; hline(2) ; 
    goodies{i} = find(zscore(std(trialsi,0,1)) < 2) ; 
end

clear stersp
for i=1:size(trigs,2) 
    ep = pop_epoch(cleaned,{trigs{i}},[-.65,2.65]) ;
    for j=1:64 ; 
        for k=1:size(ep.data,3)
            [stersp(i,j,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.data(j,:,k)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,'plotersp','off','plotitc','off',...
                'freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',NaN) ; 
        end
    end
end
stdersp = squeeze(std(mean(stersp,4),0,5)) ; 
clear goods ; for i=1:6 ; zsi = (zscore(sum(squeeze(stdersp(i,:,:)),1))) ; goods{i} = find(zsi<0.5) ; end

clear jersp ;
for i=1:size(trigs,2) 
    ep = pop_epoch(cleaned,{trigs{i}},[-.65,2.65]) ;
    for j=1:64 ; 
        [jersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.data(j,:,goods{i})),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,'plotersp','off','plotitc','off',...
            'freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',0) ; 
        
    end
end

snames = {'fovea','periphery','top hemi','bot hemi','right hemi','left hemi'} ;
nbg = find(freqs>60 & freqs<80) ; beta = find(freqs>15 & freqs<25) ; ts = find(times>0 & times<2) ; 
for i=1:6 ; subplot(2,6,i) ; topoplot(double(squeeze(mean(mean(jersp(i,:,nbg,ts),3),4))),ep.chanlocs,'maplimits',[-2,2]) ; title(snames(i)) ; end
for i=1:6 ; subplot(2,6,i+6) ; topoplot(double(squeeze(mean(mean(jersp(i,:,beta,ts),3),4))),ep.chanlocs,'maplimits',[-2,2]) ; title(snames(i)) ; end







snames = {'fovea','periphery','top hemi','bot hemi','right hemi','left hemi'} ;
for i=1:6 ; subplot(2,6,i) ; topoplot(squeeze(mean(mean(ersp(i,:,freqs>8 & freqs<15,times>0 & times<2),3),4)),EEG.chanlocs,'maplimits',[-2,2],'electrodes','off') ; title(snames{i}) ; end
for i=1:6 ; subplot(2,6,i+6) ; topoplot(squeeze(mean(mean(ersp(i,:,freqs>15 & freqs<25,times>0 & times<2),3),4)),EEG.chanlocs,'maplimits',[-2,2],'electrodes','off') ; title(snames{i}) ; end
suptitle('frequency band specific EEG retinotopy') ; 

comps = [7,9,31,34] ; 
for i=1:4 ; subplot(7,4,i) ; topoplot(squeeze(ica.icawinv(:,comps(i))),ica.chanlocs) ; end ;
icount = 4 ; 
for j=1:6
    for i=1:4 ;
        subplot(7,4,icount) ; topoplot(squeeze(ica.icawinv(:,comps(i))),ica.chanlocs) ;
        icount = icount + 1 ; 
    end 
end













    