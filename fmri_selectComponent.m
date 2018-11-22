clear all ; close all ; 
subs = {'sub_alex','sub_charest','sub_esteban','sub_fabio','sub_gab','sub_gabriella','sub_genevieve','sub_gina','sub_jeremie','sub_julie','sub_katrine','sub_lisa'...
        ,'sub_marc','sub_marie','sub_mathieu','sub_maxime','sub_mingham','sub_patricia','sub_po','sub_russell','sub_sunachakan','sub_vincent'} ;
hrf = spm_hrf(2) ; 
for sub=1:length(subs)
    cd(['c:/shared/papsaves/fmricompsubs/',subs{sub}]) ; 
    comps = load('melodic_mix') ; 
    ics = load_untouch_nii('melodic_IC.nii.gz') ; 
    icimg = ics.img ; 
    f1 = load_untouch_nii('f1.nii.gz') ; 
    icount = 1 ; 
    for i=1:245:2205
        epoched(icount,:,:) = comps(i:i+244,:) ;
        icount = icount + 1 ; 
    end
    cd trigs ; stimes = load('stimTimes1') ; stimes = stimes.stimTimes ;
    stimes = cell2mat(stimes) ; stimes = stimes(1:2:60) ; 
    trs = round(stimes/2) ; ts = zeros(1,245) ; ts(trs) = 1 ; 
    conved = conv(ts,hrf) ; conved = conved(1:245) ; 
    for i=1:size(epoched,1)
        for j=1:size(epoched,3)
            corrs(i,j) = corr2(squeeze(epoched(i,:,j)),conved) ; 
        end
    end
    mcorrs = squeeze(mean(corrs,1)) ; 
    maxinds = find(mcorrs>.29) ; allinds{sub} = maxinds ; 
    mcorrs(maxinds)
    figure,plot(mean(comps(1:245,maxinds),2),'LineWidth',2) ; suptitle(subs{sub}) ; 
    allcomps(sub,:,:) = mean(comps(:,maxinds),2) ; 
    
    stimes = dir('stimTimes*') ; 
    for st=1:length(stimes) ; 
        curr = load(stimes(st).name) ;
        curr = cell2mat(curr.stimTimes) ; 
        alltimes(st,:) = curr(1:2:60) ;
        alltrigs(st,:) = curr(2:2:60) ; 
    end
    stimcounts = [1,1,1,1,1,1] ; 
    alltimes = round(alltimes/2) ; 
    compepochs = squeeze(mean(epoched(:,:,maxinds),3)) ;
    clear stimepochs
    for i=1:size(alltimes,1) ; 
        for j=1:size(alltimes,2)
            stimij = alltrigs(i,j) ; 
            epochij = compepochs(i,(alltimes(i,j)-1):(alltimes(i,j)+8)) ; 
            stimepochs(stimij,stimcounts(stimij),:) = epochij ; 
            stimcounts(stimij) = stimcounts(stimij) + 1 ; 
        end
    end
    allstimepochs(sub,:,:,:) = stimepochs ;    
    meancomps{sub} = squeeze(mean(icimg(:,:,:,maxinds),4)) ; 
    f1.img = meancomps{sub} ; 
    save_untouch_nii(f1,'meancomp.nii.gz') ; 
end
mstim = squeeze(mean(allstimepochs,3)) ; 
stims = [1,3,2,5,6]; 
errorbar(squeeze(mean(mstim(:,stims,:),1))',squeeze(std(mstim(:,stims,:),0,1))'./sqrt(22),'LineWidth',2) ; vline([2,3],'k') ; 

for i=1:22 ; subplot(4,6,i) ; 
    errorbar(squeeze(mean(allstimepochs(i,stims,:,:),3))',squeeze(std(allstimepochs(:,stims,:,:),3))'./sqrt(45)) ; 
end
plot(mean(allcomps(:,1:245)),'LineWidth',2) ; vline(trs) ; xlabel('time(s)') ; ylabel('component strength (arb. units)') ; 
set(gca,'XTick',1:25:245,'XTickLabel',(1:25:245)*2) ; 




for i=1:length(subs)
   cd('c:/shared/allfmris/') ; cd(subs{i}) ; cd ica
   ls
   inds = allinds{i} ; 
   z150 = zeros(1,150) ; 
   z150(inds) =1 ; 
   removes = find(z150==0) ; 
   dlmwrite('removes.txt',removes) ; 
end



%%% plot the single subjects
ylabels = {'unpt-5%cont','unpt-33%cont','unpt-10%rnd','unpt-60%rnd'} ;
stypes = [2,3,5,6] ; clear tstats 
for s=1:length(stypes)
for i=1:22 ; 
    [h(s,i),p(s,i),ci,stats] = ttest(squeeze(mean(allstimepochs(i,1,:,5:6),4)),squeeze(mean(allstimepochs(i,stypes(s),:,5:6),4))) ; 
    tstats(s,i) = stats.tstat ; 
end    
end

imagesc(tstats,[-5,5]) ; xlabel('subjects(n=22)') ; set(gca,'YTick',1:4,'YTickLabel',ylabels) ; colorbar


for i=1:4 ; 
    subplot(2,2,i) ; bar(tstats(i,:)) ; 
    
    if i==1 
       xlabel('subjects(n=22), hypothesis 0/1 significant diff. (p<0.05)') ; ylabel('t-value') ;  
    end
    set(gca,'XTick',1:22,'XTickLabel',h(i,:)) ; 
    title(ylabels(i)) ; 
end

