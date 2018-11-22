cd('c:/Vision/Raw Files/qtf') ; ls  ; clear all  ; close all
sounds=dir('*vhdr') ;
for i=1:max(size(sounds)) ; 
   EEG = pop_loadbv('.',sounds(i).name) ; 
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
  % EEG = denoise_grad(EEG) ; 
   EEG = pop_resample(EEG,256) ;
   %d1 = EEG.data(1:32,:) ; EEG.data(1:32,:) = EEG.data(33:64,:) ; EEG.data(33:64,:) = d1 ; 
  % EEG = denoise_bcg(EEG) ; 
   eegs{i} = EEG ; 
   if i==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end   
end
zsums = zscore(sum(diff(merged.data,1,2).^2,2)) ;
%merged = pop_interp(merged,find(zsums>1),'spherical') ;
mergefilt = merged ; mergefilt.data = eegfiltfft(mergefilt.data,mergefilt.srate,1,128) ; 
mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,mergefilt.srate,59.5,60.5) ; 
ica = pop_runica(mergefilt,'runica') ; 

%mergeica = merged ; mergeica.icaact = icaact(merged.data,ica.icaweights*ica.icasphere,[]) ; 
%mergeica.icachansind = ica.icachansind ; mergeica.icawinv = ica.icawinv ; mergeica.icaweights = ica.icaweights ; mergeica.icasphere = ica.icasphere ; 
%%% need to get the rotation at each time point step 1 calculate single
%%% trials
comps = 1:64 ;
trigs{1} = 11:18 ; trigs{2} = 21:28 ; trigs{3} = 31:38 ; trigs{4} = 41:48 ; trigs{5} = 1:8 ; 
clear ersp 
for t1=1:length(trigs)
    for t2=1:length(trigs{t1}) ; 
        if trigs{t1}(t2) < 10
            ep = pop_epoch(ica,{['S  ',num2str(trigs{t1}(t2))]},[-3,13]) ;
        else
            ep = pop_epoch(ica,{['S ',num2str(trigs{t1}(t2))]},[-3,13]) ;
        end
        for c=1:length(comps) ; 
            for tr=1:size(ep.icaact,3)
                [ersp(t1,t2,c,tr,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(comps(c),:,tr)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                        'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'timesout',200) ; 
            end
        end
    end
end
bersp = ersp - repmat(mean(ersp(:,:,:,:,:,times<0),6),[1,1,1,1,1,200]) ; 
mbersp = squeeze(mean(mean(mean(bersp,1),2),4)) ; 
for i=1:4
    figure,
    for j=1:64 ; 
        subplot(5,13,j) ; 
        imagesc(squeeze(mean(mean(bersp(i,:,j,:,:,:),2),4)),[-4,4]) ; 
    end
    suptitle(snames{i}) ; 
end

c = [16,18] ; 
%%%% unit circle
ucirc = 1:360 ; 
% need to map each time point to an angle, or rather each angle to a time
% point because there are more angles than time points
ntimes = length(find(times>0 & times<10)) ; 
quadstarts = [225,135,45,315] ; 
rotstep = 360/ntimes ; 
for i=1:length(quadstarts) ; rotangs(i,:) = round(mod(quadstarts(i):-rotstep:quadstarts(i)-360,360)) ; end
% get the bottom left quadrant:
clear bl tr
for i=1:4
    bl{i} = find(rotangs(i,:)<270 & rotangs(i,:)>180) + length(find(times<0)) ; 
    tr{i} = find(rotangs(i,:)>0 & rotangs(i,:)<90) + length(find(times<0)) ; 
    tl{i} = find(rotangs(i,:)<180 & rotangs(i,:)>90) + length(find(times<0)) ; 
    br{i} = find(rotangs(i,:)>270 & rotangs(i,:)<360) + length(find(times<0)) ; 
end

% side note : should also compare this to static gratings in the bottom left, bottom right, etc
rbersp = squeeze(mean(mean(bersp,2),4)) ; 
for i=1:4
    meanbl(i,:,:) = squeeze(mean(rbersp(i,:,:,bl{i}),4)) ; 
    meantr(i,:,:) = squeeze(mean(rbersp(i,:,:,tr{i}),4)) ; 
    meantl(i,:,:) = squeeze(mean(rbersp(i,:,:,tl{i}),4)) ; 
    meanbr(i,:,:) = squeeze(mean(rbersp(i,:,:,br{i}),4)) ; 
end
mbl = squeeze(mean(meanbl,1)) ; 
mtr = squeeze(mean(meantr,1)) ; 
mtl = squeeze(mean(meantl,1)) ; 
mbr = squeeze(mean(meanbr,1)) ; 

plot(squeeze(mean(mbl(c,:)))) ; hold on ; plot(squeeze(mean(mbr(c,:))),'r')


angs = 1:8:360 ; 
for i=1:length(angs) ; angnames{i} = num2str(angs(i)) ; end
clear timeangles ; 
for i=1:4
    for ang=1:length(angs)-1
        timeangles{i,ang} = find(rotangs(i,:) > angs(ang) & rotangs(i,:)<angs(ang+1)) + length(find(times<0)) ; 
    end
end
clear meanangs ; 
for i=1:size(timeangles,1)
    for j=1:size(timeangles,2) ;
        meanangs(i,j,:,:) = squeeze(mean(rbersp(i,:,:,timeangles{i,j}),4)) ;    
    end
end
mmangs = squeeze(mean(meanangs,1)) ; 
for i=1:size(mmangs,1) ; 
    hold on ; 
    plot(squeeze(mean(mmangs(i,c,:),2)),'Color',[i/size(mmangs,1),0,0],'LineWidth',2) ;  
end
hline(0,'k') ; 
cc = c(1) ; 
mcangs = squeeze(mean(mmangs(:,cc,:),2)) ; 
[rho,pval] = corr(mcangs) ;
cf1 = find(freqs>10 & freqs<20) ; cf2 = find(freqs>55 & freqs<75) ; 
subplot(2,2,1) ; 
plot(squeeze(mean(mean(mmangs(:,cc,cf1),2),3)),squeeze(mean(mean(mmangs(:,cc,cf2),2),3)),'.') ; xlabel('beta') ; ylabel('gamma') ; 
title(['corr2 = ',num2str(corr2(squeeze(mean(mean(mmangs(:,cc,cf1),2),3)),squeeze(mean(mean(mmangs(:,cc,cf2),2),3))) )])
subplot(2,2,2) ; 
imagesc(freqs,angs,mcangs) ; colorbar ; xlabel('frequency(hz)') ; ylabel('angle(deg)') ; 
subplot(2,2,3) ; topoplot(ica.icawinv(:,cc),EEG.chanlocs) ;

