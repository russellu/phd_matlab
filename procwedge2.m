clear all ; close all 
cd c:/shared/retino_test/ ; ls ; 
nii = load_untouch_nii('bp_mc_Test_Russell_2015_10_22-2_WIP_EEG-fMRI_MB2_3.75mm_SENSE_7_1.nii.gz') ; nimg = nii.img ; 
cd stimtrigs ; ls 
trigs = load('scanAngle_3.mat') ; 
logfile = cell2mat(trigs.logfile) ; 
trigtimes = logfile(1:2:end) ; trigtypes = logfile(2:2:end) ; 
TR = nii.hdr.dime.pixdim(5) ; trigtimes_TR = round(trigtimes./TR) ; 
hrf = spm_hrf(TR) ; maxhrf = find(hrf==max(hrf)) ; 
maxhrfinds = maxhrf ; 
ntimes = 10/TR ; anglestep = 360/ntimes ; 
startangles = [180+45, 90+45, 0+45, 270+45] ;
for i=1:length(startangles) ; allangles(i,:) = startangles(i):-anglestep:startangles(i)-360 ; end
allangles = mod(allangles,360) ; 
%%% get the TR for every separate angle, and the baseline TR for every
%%% trial. baseline should be 12 seconds previous to the current trial
quadcounts = ones(1,4) ; clear stimvals ; 
for i=1:length(trigtimes_TR) % for all stimulus presentations
    quadind_i = floor(trigtypes(i)/10) ; 
    tr_i = trigtimes_TR(i) ; 
    for j=1:size(allangles,2)
        stimvals(:,:,:,quadcounts(quadind_i),quadind_i,j) = squeeze(mean(nimg(:,:,:,tr_i+maxhrfinds+j),4))-squeeze(mean(nimg(:,:,:,(tr_i+maxhrfinds+j)-round(8/TR)),4)) ; 
        sbase(:,:,:,quadcounts(quadind_i),quadind_i,j) = squeeze(mean(nimg(:,:,:,(tr_i+maxhrfinds+j)-round(8/TR)),4)) ; 
        sstim(:,:,:,quadcounts(quadind_i),quadind_i,j) = squeeze(mean(nimg(:,:,:,tr_i+maxhrfinds+j),4)) ; 
    end
    quadcounts(quadind_i) = quadcounts(quadind_i) + 1 ; 
end

% do the t-test for each angle
tbrains = zeros(size(squeeze(mean(stimvals,4)))) ; 
for i=1:size(stimvals,5) ; 
    for j=1:size(stimvals,6)
        sbaseij = squeeze(sbase(:,:,:,:,i,j)) ; 
        sstimij = squeeze(sstim(:,:,:,:,i,j)) ; 
        resbase = reshape(sbaseij,[size(sbase,1)*size(sbase,2)*size(sbase,3),size(sbase,4)]) ; 
        resstim = reshape(sstimij,[size(sstimij,1)*size(sstimij,2)*size(sstimij,3),size(sstimij,4)]) ; 
        [h,p,ci,tstat] = ttest(resbase',resstim') ; 
        tstats = tstat.tstat ; 
        tbrainij = reshape(tstats,[size(sbaseij,1),size(sbaseij,2),size(sbaseij,3)]) ; 
        tbrains(:,:,:,i,j) = tbrainij ; 
    end
end
%mstimvals = squeeze(mean(stimvals,4))./squeeze(std(stimvals,0,4)) ; 
mstimvals = tbrains ; 
% need to start the first stimulus a bit later...
% assign an angle to every index, and then sort them?
% you need to assign a linear index to every i,j in the matrix
% linear index from 1-44 or whatever
for i=1:size(allangles,1)
    for j=1:size(allangles,2)
        inds(i,j) = sub2ind(size(allangles),i,j) ; 
    end
end
% do a t-test for each of the n stimulus presentations.
resangles = reshape(allangles,[1,numel(allangles)]) ; 
resinds = reshape(inds,[1,numel(allangles)]) ; 
[sv,si] = sort(resangles,'descend') ; 
sortinds = resinds(si) ; % sortinds are the indices for each angle sorted according to the degree of the angle
[cx,cy] = ind2sub(size(allangles),sortinds) ; 
angleinds = zeros(size(allangles)) ; 
for i=1:length(cx) ;
    angleinds(cx(i),cy(i)) = i ; 
end

orderstimvals = zeros(size(mstimvals,1),size(mstimvals,2),size(mstimvals,3),numel(angleinds)) ; 
mstimvals(isnan(mstimvals)) = 0 ; mstimvals(isinf(mstimvals)) = 0 ; 
restim = reshape(mstimvals,[1,numel(mstimvals)]) ; restim(abs(zscore(restim))>5) = 0 ; mstimvals = reshape(restim,size(mstimvals)) ; 
for i=1:size(mstimvals,4) 
    for j=1:size(mstimvals,5) 
        orderstimvals(:,:,:,angleinds(i,j)) = mstimvals(:,:,:,i,j) ; 
    end
end

%for i=1:36 ; subplot(6,6,i) ; imagesc(squeeze(mean(orderstimvals(:,:,8:12,i),3))) ; end

smoothvals = zeros(size(orderstimvals)) ; 
for i=1:size(orderstimvals,1)
    for j=1:size(orderstimvals,2)
        for k=1:size(orderstimvals,3)
            smoothvals(i,j,k,:) = smooth(squeeze(orderstimvals(i,j,k,:))) ; 
        end
    end
end
for i=1:size(smoothvals,4)  
    smoothvals(:,:,:,i) = medfilt3(squeeze(smoothvals(:,:,:,i))); 
end
save_nii(make_nii(smoothvals),'smoothvals.nii.gz') ; 

stdmask = squeeze(std(smoothvals,0,4)) ; 
stdmask = imdilate(stdmask>.8,strel(ones(7,7,7))) ; 
postmask = zeros(size(stdmask)) ; 
postmask(:,48:end,:) = 1 ; 

clear cx cy cz ; 
maskret = smoothvals.*repmat(postmask,[1,1,1,size(smoothvals,4)]).*repmat(stdmask,[1,1,1,size(smoothvals,4)]) ; 
centmask = zeros(size(smoothvals)) ; 
for i=1:size(maskret,4)
    voli = squeeze(maskret(:,:,:,i)) ; 
    newvoli = zeros(size(voli)) ; 
    greaters = voli(voli>0) ; greaterinds = find(voli>0) ; 
    [sv,si] = sort(greaters,'descend') ;
    topinds = greaterinds(si(1:25)) ; 
    newvoli(topinds) = 1 ; %newvoli = medfilt3(newvoli) ; 
    maskret(:,:,:,i) = newvoli ;    
    [cx(i),cy(i),cz(i)] = centmass3(newvoli) ;
end
subplot(1,3,1) ; plot((cx)) ; subplot(1,3,2) ; plot((cy)) ; subplot(1,3,3) ; plot((cz)) ;  
%binsmooth = (smoothvals>1.5).*repmat(postmask,[1,1,1,44]) ; 

save_nii(make_nii(maskret),'maskret.nii.gz') ; 

%{
%%%% do the EEG part
cd('C:\shared\badger\russ\russ') ; ls ; clear all ; close all ; 
%filenames = {'Retino_angles_2.vhdr','Retino_angles_3.vhdr','Retino_angles_4.vhdr','Retino_angles_5.vhdr','Retino_angles_6.vhdr','Retino_polar_1.vhdr','rest.vhdr'} ; 
filenames=dir('*vhdr') ; 
for fname = 1:length(filenames)
    EEG = pop_loadbv('.',filenames(fname).name) ; 
    EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
    EEG = denoise_grad2(EEG) ; 
    EEG = denoise_bcg(EEG) ; 
    eegs{fname} = EEG ; 
    if fname==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end   
end
mergefilt = merged ; mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,mergefilt.srate,0,0.5) ; 
mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,mergefilt.srate,59.5,60.5) ; 
ica = pop_runica(mergefilt,'runica') ;

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
        imagesc(squeeze(mean(mean(bersp(i,:,j,:,:,:),2),4)),[-2,2]) ; 
    end
end

mbersp = squeeze(mean(mean(mean(bersp(1:4,:,:,:,:,:),1),2),4)) ; 
for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mbersp(i,:,:)),[-2,2]) ; title(i) ; end

c = [21,22,34] ; 
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
cc = c(3) ; 
mcangs = squeeze(mean(mmangs(:,cc,:),2)) ; 
[rho,pval] = corr(mcangs) ;
cf1 = find(freqs>10 & freqs<20) ; cf2 = find(freqs>55 & freqs<75) ; 
subplot(2,2,1) ; 
plot(squeeze(mean(mean(mmangs(:,cc,cf1),2),3)),squeeze(mean(mean(mmangs(:,cc,cf2),2),3)),'.') ; xlabel('beta') ; ylabel('gamma') ; 
title(['corr2 = ',num2str(corr2(squeeze(mean(mean(mmangs(:,cc,cf1),2),3)),squeeze(mean(mean(mmangs(:,cc,cf2),2),3))) )])
subplot(2,2,2) ; 
imagesc(freqs,angs,mcangs) ; colorbar ; xlabel('frequency(hz)') ; ylabel('angle(deg)') ; 
subplot(2,2,3) ; topoplot(ica.icawinv(:,cc),EEG.chanlocs) ;

%%%%% correlate the power at that angle at that electrode with the distance
%%%%% from the center of mass of the activation cluster.
%%% still need a way to calculate power on a trial by trial basis...that doesn't use baseline correction or taking the log...after ICA cleaning,
%%% so you can have power values in each electrode
%}
