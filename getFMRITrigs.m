cd('C:\shared\badger\russ\russ') ; 
clear all ; close all ; 
mongs = dir('*vhdr') ; 
for i=1:length(mongs) ; 
   EEG = pop_loadbv('.',mongs(i).name) ;  
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
   EEG = denoise_grad2(EEG) ; 
   EEG = denoise_bcg(EEG) ;    
   eegs{i} = EEG ; 
   if i==1 ; merged = EEG  ;else merged = pop_mergeset(EEG,merged) ; end     
end
mergefilt = merged ; mergefilt.data = eegfiltfft(mergefilt.data,mergefilt.srate,1,128) ; 
mergefilt.data = mergefilt.data - eegfiltfft(mergefilt.data,mergefilt.srate,59.5,60.5) ; 
ica = pop_runica(mergefilt,'runica') ; 

goodcs = [10,24,26,31] ; allcs = zeros(1,64) ; allcs(goodcs) = 1 ;
ica2 = pop_subcomp(ica,find(allcs==0)) ; 
comps = 1:64 ;
trigs{1} = 11:18 ; trigs{2} = 21:28 ; trigs{3} = 31:38 ; trigs{4} = 41:48 ; trigs{5} = 1:8 ; 
clear ersp 
for t1=1:length(trigs)
    for t2=1:length(trigs{t1}) ; 
        if trigs{t1}(t2) < 10
            ep = pop_epoch(ica2,{['S  ',num2str(trigs{t1}(t2))]},[-3,13]) ;
        else
            ep = pop_epoch(ica2,{['S ',num2str(trigs{t1}(t2))]},[-3,13]) ;
        end
        for c=1:length(comps) ; 
            for tr=1:size(ep.icaact,3)
                [ersp(t1,t2,c,tr,:,:),itc,powbase,times,freqs,~,~] = nonlogtimef(squeeze(ep.data(comps(c),:,tr)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                        'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'timesout',200) ; 
            end
        end
    end
end
bersp = ersp - repmat(mean(ersp(:,:,:,:,:,times<0),6),[1,1,1,1,1,200]) ; 
mbersp = squeeze(mean(mean(mean(bersp,1),2),4)) ; 

elabs = {ica.chanlocs.labels} ;
for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mersp(i,:,:)),[-.2,.2]) ; title(elabs{i}) ; end
%c = [24,26] ; 
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
%{
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
%}

for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mmangs(:,i,:)),[-.1,.1]) ; title(elabs{i}) ; end

%%% get the triggers from the EEG files for locating the BOLD times
nonstims = {'Sync On','boundary','R128','S 98','S 99'} ; 
cd('C:\shared\badger\russ\topup') ;
% get the triggers first:
for eg=2:6
    template = load_untouch_nii('template.nii.gz') ; 
    stimtriginds = [] ; stimtrigs = [] ; 
    eegi = eegs{eg} ;
    evts = {eegi.urevent.type} ;
    r128s = find(strcmp(evts,'R128')) ; 
    lats = {eegi.urevent.latency} ; 
    % need to find the volume closest to each trigger.
    fmri = load_untouch_nii(['c_reg_mc_Test_Russell_2015_10_15_WIP_EEG-fMRI_MB2_3mm_SENSE_1',num2str(eg),'_1.nii.gz']) ; 
    for e=1:length(evts)  ;
        index = cell2mat(strfind(nonstims,evts{e})) ; 
        if sum(index) == 0 
            stimtrigs(length(stimtrigs)+1) = str2num(evts{e}(3:4)) ; 
            stimtriginds(length(stimtriginds)+1) = e ; 
        end
    end
    stimlats = cell2mat(lats(stimtriginds)) ; 
    r128lats = cell2mat(lats(r128s)) ; 
    % find the closest BOLD volume to each stimulus trigger
    for s=1:length(stimlats)
        latdiffs = abs(stimlats(s) - r128lats) ; 
        mindiffinds(s) = find(latdiffs == min(latdiffs)) ; 
    end
       
    trigtypes = stimtrigs ; 
    TR = fmri.hdr.dime.pixdim(5) ; 
    trigtimes_TR = mindiffinds ; 
    
    nimg = fmri.img ; 
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
    anglevals = zeros(1,numel(allangles)) ; 
    for i=1:size(mstimvals,4) 
        for j=1:size(mstimvals,5) 
            orderstimvals(:,:,:,angleinds(i,j)) = mstimvals(:,:,:,i,j) ; 
            anglevals(angleinds(i,j)) = allangles(i,j) ; 
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
    template.img = smoothvals ; 
    save_untouch_nii(template,['smoothvals_',num2str(eg),'.nii.gz']) ;     
end

% need a correlation metric...distance vs BOLD-EEG relationship?
% tuning curve in each voxel vs distance
% correlation coefficient of each voxel's tuning curve, vs distance?
% get the mask
smoothie = load_untouch_nii('allsmooth_t1.nii.gz') ; smoothim = smoothie.img ; 
smoothstd = (squeeze(std(smoothim,0,4))) ;
t1 = load_untouch_nii('../t1.nii') ; 
[kc,km] = kmeans(uint8(mat2gray(reshape(smoothstd,[1,numel(smoothstd)]))*255),3) ; 
kimg = reshape(km,size(smoothstd)) ; 
locs = load_untouch_nii('../t1_locs.nii.gz') ; 
elecorder = load('../elecorder') ; elecorder = elecorder.elecorder ; 
elabs = {ica.chanlocs.labels} ; 
clear cx cy cz ; 
for i=1:length(elecorder)
    [cx(i),cy(i),cz(i)] = centmass3(locs.img==i) ; 
end
visinds = find(kimg==3) ;
[vx,vy,vz] = ind2sub(size(kimg),visinds) ; 
for i=1:length(cx) ; 
    vxdiffs(i,:) = cx(i)-vx ; 
    vydiffs(i,:) = cy(i)-vy ; 
    vzdiffs(i,:) = cz(i)-vz ; 
end
% get the center of mass of each cluster
for i=1:size(smoothim,4)
    smoothi = squeeze(smoothim(:,:,:,i)) ; 
    visi = smoothi(visinds) ; 
    gthresh = visi>6 ; 
    centvis = zeros(size(smoothi)) ; 
    centvis(visinds(gthresh==1)) = 1 ; 
    %t1.img = centvis ; save_untouch_nii(t1,['centvis',num2str(i),'.nii.gz']) ; 
    cents(:,:,i) = squeeze(sum(centvis,1)) ; 
    [bcx(i),bcy(i),bcz(i)] = centmass3(centvis) ; 
end

for i=1:length(cx)
    for j=1:length(bcx)
        elecdists(i,j) = sqrt((cx(i)-bcx(j)).^2 + (cy(i)-bcy(j)).^2 + (cz(i)-bcz(j)).^2)  ;                
    end 
end

bangles = squeeze(mean(mmangs(:,:,freqs<25 & freqs>12),3)) ; 
% for each electrode, a scatter plot of distance vs power
% start with poz.
clear revangles 
for i=1:64 ; revangles(:,i) = abs(smooth(flipud(imresize(bangles(:,i),[36,1])))) ; end
clear vals ; 
for i=1:length(elecorder)
   eegind = find(strcmpi(elecorder{i},elabs)) ;
   if ~isempty(eegind)
       corrs(i) = corr2(revangles(:,eegind),elecdists(i,:)') ; 
       vals(i,:,1) = revangles(:,eegind) ;
       vals(i,:,2) = elecdists(i,:)' ;
   end
end

for i=1:length(elabs) ; if i~=32 ; orderinds(i) = find(strcmpi(elabs{i},elecorder)) ; end ; end
locinds = orderinds(orderinds~=0) ; labinds = [1:31,33:64] ; 
clear corrs vals ; 
for i=1:size(revangles,1)
    corrs(i) = corr2(elecdists(locinds,i),revangles(i,labinds)') ; 
    vals(i,:,1) = elecdists(locinds,i) ;
    vals(i,:,2) = revangles(i,labinds)' ; 
end





