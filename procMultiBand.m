clear all ; close all ; 
cd ('c:/shared/Test MB-SENSE 6') ; ls 
tr6 = load_nii('Test_MB-SENSE_6_WIP_fMRI_S2_MB3_SENSE_8_1.nii') ;
%tr6 = load_nii('Test_MB-SENSE_6_WIP_fMRI_S2_MB1_SENSE_6_1.nii') ;
trim = tr6.img ; 
st1 = load('stimTimes1.mat') ; st2 = load('stimTimes2.mat') ; 
st1 = st1.stimTimes ; st2 = st2.stimTimes ; 
tr = .8 ; 
clear times trs labs
for i=1:size(st1,2)
    times(i) = st2{i}(1) ; 
    trs(i) = round(times(i)/tr) ; 
    labs(i) = st2{i}(2) ; 
end
block = zeros(1,size(trim,4)) ; 
hrf = spm_hrf(tr) ; 
block(trs) = 1 ; 
convblock  = conv(block,hrf) ; 
convblock = convblock(1:size(trim,4)) ; 
dlmwrite('convblock2',convblock') ;
% do stuff in afni (get the corrmask)
corr = load_nii('conv_corrs.nii.gz') ; corrimg = corr.img ; 
postmask = (medfilt3(corr.img>.10)) ; 
stims = [1,2,3,4,5,6] ; 
stimnames = {'unpterb.','5%cont.','33%cont','plaid','10%rnd','60%rnd'} ; 
postinds = find(postmask==1) ; 
for i=1:max(size(postinds))
    [i1,i2,i3] = ind2sub(size(postmask),postinds(i)) ; 
    voxels(i,:) = squeeze(trim(i1,i2,i3,:)) ; 
    corrs(i) = corrimg(i1,i2,i3) ; 
end
[sortv,sorti] = sort(corrs,'descend') ; 
sortvoxels(1:size(voxels,1),:) = voxels(sorti,:) ; 
subplot(3,1,1) ; plot(squeeze(mean(sortvoxels(1:50,:))))
msortvoxels = squeeze(mean(sortvoxels(1:100,:),1)) ; 

%%% define the time periods

basel = round(4/tr) ; 
taskl = round(15/tr) ; 
clear stimepochs normepochs
for st=1:max(size(stims))
sinds = find(labs==stims(st)) ; 
sindts = trs(sinds) ; 
for i=1:size(sindts,2)
    stimepochs(st,i,:) = msortvoxels(sindts(i)-basel:sindts(i)+taskl) ; 
    normepochs(st,i,:) = (stimepochs(st,i,:)-squeeze(mean(stimepochs(st,i,1:basel),3)))./squeeze(mean(stimepochs(st,i,1:basel),3)) ; 
end
end
timess = -basel*tr:tr:taskl*tr ; 
subplot(2,2,1) ; 
errorbar(squeeze(mean(normepochs(1,:,:),2)),squeeze(std(normepochs(1,:,:),0,2))./sqrt(5),'r','LineWidth',2) ; hold on ; 
errorbar(squeeze(mean(normepochs(2,:,:),2)),squeeze(std(normepochs(1,:,:),0,2))./sqrt(5),'g','LineWidth',2) ;
errorbar(squeeze(mean(normepochs(3,:,:),2)),squeeze(std(normepochs(1,:,:),0,2))./sqrt(5),'b','LineWidth',2) ;
title('3 contrast lvls') ; xlabel('time(s)') ; ylabel('(task-rest)/rest') ; set(gca,'XTick',1:2:size(normepochs,3),'XTickLabel',timess(1:2:end)) ; hline(0,'k') ; ylim([-.01,.045]) ; legend(stimnames{[1,2,3]}) ; xlim([1,size(stimepochs,3)]) ;
vline(find(timess==0)) ; 

subplot(2,2,2) ;
errorbar(squeeze(mean(normepochs(1,:,:),2)),squeeze(std(normepochs(1,:,:),0,2))./sqrt(5),'r','LineWidth',2) ; hold on ; 
errorbar(squeeze(mean(normepochs(5,:,:),2)),squeeze(std(normepochs(1,:,:),0,2))./sqrt(5),'c','LineWidth',2) ;
errorbar(squeeze(mean(normepochs(6,:,:),2)),squeeze(std(normepochs(1,:,:),0,2))./sqrt(5),'m','LineWidth',2) ;
title('3 randomization lvls') ; xlabel('time(s)') ; ylabel('(task-rest)/rest') ; set(gca,'XTick',1:2:size(normepochs,3),'XTickLabel',timess(1:2:end)) ; hline(0,'k') ; ylim([-.01,.045]) ; legend(stimnames{[1,5,6]}) ;  xlim([1,size(stimepochs,3)]) ;
vline(find(timess==0)) ; 

normepochs = double(normepochs) ; 
%%%% do some stats
for i=1:size(normepochs,3) ;
    [p(i),atab,stats] = anova1([squeeze(normepochs(1,:,i));squeeze(normepochs(2,:,i))]',[],'off') ; 
    f(i) = atab{2,5} ; 
    [p2(i),atab,stats] = anova1([squeeze(normepochs(1,:,i));squeeze(normepochs(6,:,i))]',[],'off') ; 
    f2(i) = atab{2,5} ; 
end
figure,
subplot(2,2,1) ; plot(p) ; 
title('100% contrast > 5% contrast') ; xlabel('time(s)') ; ylabel('p-value') ; set(gca,'XTick',1:2:size(normepochs,3),'XTickLabel',timess(1:2:end)) ;
ylim([0,1]) ; xlim([1,size(p,2)]) ; vline(find(timess==0)) ; 
subplot(2,2,2) ; plot(p2) ; 
title('60% rnd > 0% rnd') ; xlabel('time(s)') ; ylabel('p-value') ; set(gca,'XTick',1:2:size(normepochs,3),'XTickLabel',timess(1:2:end)) ;
ylim([0,1]) ; xlim([1,size(p,2)]) ; vline(find(timess==0)) ; 







%%%%%% epoch and average the data or do a 2d correlation in every voxel
%%%%%% (stimulus types)




