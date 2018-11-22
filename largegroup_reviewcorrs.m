clear all ;close all ; 
subs = {'a_alex','a_charest','a_esteban','a_fabio','a_gab','a_gabriella','a_genevieve','a_gina','a_guillaume','a_jeremie','a_julie','a_katrine','a_lisa','a_marc',...
    'a_marie','a_mathieu','a_maxime','a_mingham','a_patricia','a_po','a_russell','a_sunachakan','a_tah','a_vincent'} ; 

eeg = pop_loadset('e:/nimg_pool/a_alex/cleanfilt.set') ; 


for sb=1:length(subs)
    
    % DISTANCE
    cd(['e:/nimg_pool/',subs{sb}]) ;  
    
    corrs = load_untouch_nii('cleancorrs_fs.nii.gz'); 
    coords = load('coords.mat'); coords = coords.coords; 
    size(coords); 
    subcoords(sb,:,:) = coords; 
    
    % create ROI
    corrs.img(:,:,100:end) = 0; 
    zcorrs = zeros(size(corrs.img)); 
    corrs.img(isnan(corrs.img)) = 0; 
    [sv,si] = sort(corrs.img(:),'descend'); 
    zcorrs(si(1:15000)) = 1; 
    
    % get distance from ROI
    [cx,cy,cz] = centmass3(zcorrs); 
    sx = size(corrs.img,1); sy = size(corrs.img,2); sz = size(corrs.img,3); 
   % [xg,yg,zg] = ndgrid(-sx/2:sx/2-1,-sy/2:sy/2-1,-sz/2:sz/2-1);
    [xg,yg,zg] = ndgrid(-cx:sx-cx-1,-cy:sy-cy-1,-cz:sz-cz-1);
    sphere = sqrt(xg.^2 + yg.^2 + zg.^2) < 5; 
    
    corrsphere = corrs ; corrsphere.img = double(sphere) ; save_untouch_nii(corrsphere,'corrsphere.nii.gz'); 
    
    coord_diffs = [coords(:,1) - cx, coords(:,2) - cy, coords(:,3) - cz];
    dists = sqrt(sum(coord_diffs.^2,2)); 
    subdists(sb,:) = dists; 

    
    brain = load_untouch_nii('fs_t1.nii.gz'); 
    %zmask  = double(zcorrs> 0 | corrsphere.img>0); 
    %maskcorrs = zmask; maskcorrs(sphere>0) = -1; 
    %subplottight(2,12,sb);plotoverlayIntensity2D(squeeze(brain.img(:,cy,:)),squeeze(zmask(:,cy,:)),squeeze(maskcorrs(:,cy,:)),90); 

    leftnorms = load_untouch_nii('left_norms.nii.gz') ; 
    rightnorms = load_untouch_nii('right_norms.nii.gz') ; 
    bothnorms = leftnorms.img + rightnorms.img ; 
    
    norm_mask = sum(abs(bothnorms),4)>0; 
    norm_corrs = corrs.img.*norm_mask;
    norm_corrs(isnan(norm_corrs)) = 0; 
    [sv,si] = sort(norm_corrs(:),'descend'); 
    
    normls = 100:250:5000 ; 
    for nl=1:length(normls)  
        [nx,ny,nz] = ind2sub(size(zcorrs),si(1:normls(nl))) ; 
        [bx,by,bz] = ind2sub(size(zcorrs),si(1:2500)); 
        allnorms = zeros(length(nx),3) ; 
        allbnorms = zeros(length(bx),3); 
        for i=1:length(nx)
            allnorms(i,:) = squeeze(bothnorms(nx(i),ny(i),nz(i),:)) ; 
        end
        for i=1:length(bx)
            allbnorms(i,:) = squeeze(bothnorms(bx(i),by(i),bz(i),:)) ; 
        end
        binots(sb) = 1 - norm(sum(allbnorms,1))./length(bx) ; 
        inots(sb,nl) = 1 - norm(sum(allnorms,1))/length(nx); 
        
        
        subnorms(sb,:) = (sum(allbnorms,1))./length(bx); 
    end
  
end


cd E:\nimg_pool\saved
times = load('times'); times = times.times ; freqs = load('freqs') ; freqs = freqs.freqs; 
a_elabs = load('a_elabs'); a_elabs = a_elabs.a_elabs; 
elecorder = load('elecorder'); elecorder = elecorder.elecorder; 
all_postchans = {'O1','OZ','OZ','PO3','POZ','PO4','P8','P6','P4','P2','PZ','P1','P3','P5','P7'}; 
elecns = [23,56,24,57,25,58,26,59,27,60,61,62,63,64,29,30,31];
for i=1:length(all_postchans)
   a_chaninds(i) = find(strcmpi(all_postchans{i},a_elabs)); 
   elecorder_chaninds(i) = find(strcmpi(all_postchans{i},elecorder)); 
end
for i=1:length(a_elabs)
    ind = find(strcmpi(a_elabs{i},elecorder));
    if ~isempty(ind)
    all_chaninds(i) = find(strcmpi(a_elabs{i},elecorder));
    end
end
elecdists = zeros(24,64) ; for i=1:length(all_chaninds) ; if all_chaninds(i) ~= 0 ; elecdists(:,i) = subdists(:,all_chaninds(i)); end ; end
plotchans = zeros(1,64) ; plotchans([46,22,17,41]) = 1; goods = find(plotchans==0); 



allsub_allstim_alpha = load('allsub_allstim_alpha'); alpha = allsub_allstim_alpha.allsub_allstim_alpha; 
allsub_allstim_gamma = load('allsub_allstim_gamma'); gamma = allsub_allstim_gamma.allsub_allstim_gamma; 
allstim_ersp = load('allstim_ersp'); allstim_ersp = allstim_ersp.allstim_ersp; 
mersp = squeeze(mean(allstim_ersp(:,1:6,:,:,:),2)); 

mgamma = squeeze(mean(gamma(:,:,a_chaninds),3)); 
malpha = squeeze(mean(alpha(:,:,a_chaninds),3)); 
[eleccorrs,elecps] = corr(mean(mgamma,2),elecdists);

% figure 1
mdists = squeeze(mean(subdists(:,elecorder_chaninds),2)); 
continds = [1,3,2]; rndinds = [5,4,6]; contlabs = {'100% contrast','5% contrast','33% contrast'}; rndlabs = {'10% random','plaid','60% random'}; 
stimlabs = {contlabs{1},contlabs{2},contlabs{3},rndlabs{1},rndlabs{2},rndlabs{3}}; 
contcolors = {[0.99,0,0],[0.5,0,0.25],[0.25,0,0.5]}; rndcolors = {[0,0.0,0.99],[0,0.25,0.5],[0,0.5,0.25]}; 
stimcolors = {contcolors{1},contcolors{2},contcolors{3},rndcolors{1},rndcolors{2},rndcolors{3}};


fig1=figure,
subplot(4,8,1); 
axis([0,120,-4,1.5]); hold on ; redbox = [1,0.85,0.85]; bluebox = [0.85,0.85,1];
rectangle('Position',[8,-2.4,17,6],'FaceColor',bluebox,'EdgeColor',bluebox)
rectangle('Position',[40,-2.4,60,6],'FaceColor',redbox,'EdgeColor',redbox)
for i=1:3 
    shadedErrorBar(1:2:120,squeeze(mean(mean(mersp(:,continds(i),:,35:80),1),4)),squeeze(std(mean(mersp(:,continds(i),:,35:80),4),0,1))/sqrt(24),{'Color',contcolors{i}}); hold on ; hline(0,'k'); xlabel('frequency(hz)'); ylabel('dB'); 
end
ylim([-2.5,.8]); 
title('contrast');
subplot(4,8,2); 
axis([0,120,-4,1.5]); hold on ; redbox = [1,0.85,0.85]; bluebox = [0.85,0.85,1];
rectangle('Position',[8,-2.4,17,6],'FaceColor',bluebox,'EdgeColor',bluebox)
rectangle('Position',[40,-2.4,60,6],'FaceColor',redbox,'EdgeColor',redbox)
for i=1:3 
    shadedErrorBar(1:2:120,squeeze(mean(mean(mersp(:,rndinds(i),:,35:80),1),4)),squeeze(std(mean(mersp(:,rndinds(i),:,35:80),4),0,1))/sqrt(24),{'Color',rndcolors{i}}); hold on ; hline(0,'k'); xlabel('frequency(hz)'); ylabel('dB');
end
ylim([-2.5,.8]); 
title('randomization, plaid'); 
%'emarker2',{elecns,'.','r'}
subplot(4,8,3) ; topoplot(squeeze(mean(mean(gamma,1),2)),eeg.chanlocs,'electrodes','off','style','both','emarker2',{elecns,'.','r'}); title('gamma (40-100Hz)'); h = colorbar; title(h,'\muV'); colormap parula %h=colorbar ; set(get(h,'label'),'string','\muV');
subplot(4,8,4) ; topoplot(squeeze(mean(mean(alpha,1),2)),eeg.chanlocs,'electrodes','off','style','both','emarker2',{elecns,'.','r'}); title('alpha/beta (8-25Hz)'); h = colorbar; title(h,'\muV'); colormap parula
subplot(4,8,5) ; imagesc(times,freqs,squeeze(mean(mersp(:,1,:,:),1)),[-3,1.5]) ; axis xy  ; xlabel('time(s)'); ylabel('frequency(hz)'); colorbar ; vline([0,2],'k');  h = colorbar; title(h,'dB'); 
subplot(4,8,9) ; plot(1,'Color',contcolors{1},'LineWidth',5); hold on; plot(2,'Color',contcolors{2},'LineWidth',5) ; plot(3,'Color',contcolors{3},'LineWidth',5); legend(contlabs);
subplot(4,8,10) ; plot(1,'Color',rndcolors{1},'LineWidth',5); hold on; plot(2,'Color',rndcolors{2},'LineWidth',5) ; plot(3,'Color',rndcolors{3},'LineWidth',5); legend(rndlabs); 
figure,
for i=1:24 ; subplottight(2,12,i) ; imagesc(squeeze(mean(mersp(i,[1,5],:,:),2)),[-3,3]) ; axis xy ;  set(gca,'XTick',[],'YTick',[]); text(10,58,['subject ',num2str(i)]); end
set(fig1,'Units','normalized');
set(fig1,'Position',[0 0 1.2 .8]);

%%%% FIGURE 2

[freq_corrs,freq_ps] = corr(squeeze(mean(mean(mersp(:,:,:,times>0 & times<2),2),4))); 

fig2=figure,

mersp_alpha = squeeze(mean(mean(mersp(:,:,4:12,times>0.5 & times<2),3),4)); 
mersp_gamma = squeeze(mean(mean(mersp(:,:,20:50,times>0.5 & times<2),3),4)); 

subplot(4,8,1); 
for i=2:6
    plot(malpha(:,1),malpha(:,i),'d','Color',stimcolors{i},'LineWidth',3); hold on ; lsline ; xlabel('\muV alpha/beta (8-25 Hz)'); ylabel('\muV alpha (8-25 Hz)'); 
    [rhos(i),ps(i)] = corr(malpha(:,1),malpha(:,i)); 
end
title(['median ',format_rho(median(rhos)),' median ',(format_p(median(ps)))]);

subplot(4,8,2); 
for i=2:6
    plot(mgamma(:,1),mgamma(:,i),'d','Color',stimcolors{i},'LineWidth',3); hold on ; lsline ;  xlabel('\muV gamma (40-100 Hz)'); ylabel('\muV gamma (40-100-25 Hz)'); 
    [rhos(i),ps(i)] = corr(mgamma(:,1),mgamma(:,i)); 
end
title(['median ',format_rho(median(rhos)),' median ',(format_p(median(ps)))]);

subplot(4,8,3); 
for i=1:6
    plot(malpha(:,i),mgamma(:,i),'d','Color',stimcolors{i},'LineWidth',3); hold on ; lsline ;  xlabel('\muV alpha/beta (8-25 Hz)'); ylabel('\muV gamma (40-100 Hz)'); 
    [rhos(i),ps(i)] = corr(malpha(:,i),mgamma(:,i)); 
end
title(['median ',format_rho(median(rhos)),' median ',(format_p(median(ps)))]);

subplot(4,8,4);
%plot(squeeze(mean(mgamma,1)),mean(malpha,1),'o'); lsline ; hold on ; 
errorbarxy(mean(malpha,1),mean(mgamma,1),std(malpha,0,1)/sqrt(24),std(mgamma,0,1)/sqrt(24),{'.k','k','k'}); lsline ; hold on ; 
for i=1:6 ; plot(squeeze(mean(malpha(:,i),1)),squeeze(mean(mgamma(:,i),1)),'o','Color',stimcolors{i},'LineWidth',3); end
[rho,p] = corr(mean(mgamma,1)',mean(malpha,1)'); 
title([format_rho((rho)),' ',(format_p((p)))]); xlabel('\muV alpha/beta (8-25Hz)'); ylabel('\muV gamma (40-100Hz)')

subplot(4,8,5) ; for i=1:6; plot(i,'Color',stimcolors{i},'LineWidth',3) ; hold on ; end ; legend(stimlabs); 


for i=1:6 ; subplot(4,8,8+i); 
    imagesc(times,freqs,squeeze(mersp(21,i,:,:)),[-3,3]); axis xy ; 
        if i==1 ;xlabel('time(s)') ;ylabel('frequency(hz)'); end

    subplot(4,8,16+i); 
    imagesc(times,freqs,squeeze(mersp(14,i,:,:)),[-3,3]); axis xy ; 
end
subplot(4,8,16) ; bar(squeeze(mean(mean(mersp(21,:,20:50,times>0 & times<2),3),4)),'r'); hold on ; bar(squeeze(mean(mean(mersp(21,:,5:12,times>0 & times<2),3),4)),'b'); ylim([-3,1.5]); xlim([0.5,6.5])
legend({'gamma','alpha/beta'}) ; xlabel('stimulus type') ; ylabel('dB'); 
subplot(4,8,24) ; bar(squeeze(mean(mean(mersp(14,:,20:50,times>0 & times<2),3),4)),'r'); hold on ; bar(squeeze(mean(mean(mersp(14,:,5:12,times>0 & times<2),3),4)),'b'); ylim([-3,1.5]); xlim([0.5,6.5])

subplot(4,8,26) ; topoplot(squeeze(mean(gamma(21,:,:),2)),eeg.chanlocs,'maplimits',[-.1,0.1],'style','map') ; colormap parula ; h=colorbar ; title(h,'\muV gamma'); 
subplot(4,8,27) ; topoplot(squeeze(mean(alpha(21,:,:),2)),eeg.chanlocs,'maplimits',[-1,1],'style','map') ; colormap parula ; h=colorbar ; title(h,'\muV alpha/beta'); 
subplot(4,8,28) ; topoplot(squeeze(mean(gamma(14,:,:),2)),eeg.chanlocs,'maplimits',[-.1,0.1],'style','map') ; colormap parula ; h=colorbar ; title(h,'\muV gamma'); 
subplot(4,8,29) ; topoplot(squeeze(mean(alpha(14,:,:),2)),eeg.chanlocs,'maplimits',[-1,1],'style','map') ; colormap parula ; h=colorbar ; title(h,'\mu alpha/beta'); 

subplot(4,8,30) ; imagesc([-3,3]); h=colorbar; title(h,'dB'); 

set(fig2,'Units','normalized');
set(fig2,'Position',[0 0 1.2 .8]);

% FIGURE 3
fig3=figure,
subplot(4,8,1);
[rho,p] = corr(mean(mgamma,2),mdists); 
plot(mdists,squeeze(mean(mgamma(:,:),2)),'kd','LineWidth',2); lsline ; xlabel('distance (mm)'); ylabel('\muV gamma'); 
title([format_rho(rho),' ',format_p(p)]);
subplot(4,8,2); 
[rho,p] = corr(mean(malpha,2),mdists); 
plot(squeeze(mean(malpha(:,:),2)),mdists,'kd','LineWidth',2); lsline ; xlabel('distance (mm)'); ylabel('\muV alpha/beta '); 
title([format_rho(rho),' ',format_p(p)]);

subplot(4,8,3);
topoplot(eleccorrs,eeg.chanlocs,'maplimits',[-.6,0],'plotchans',goods,'style','map') ; colormap parula  ;h=colorbar ;title(h,'rho'); title('gamma vs distance'); 
subplot(4,8,4); 
topoplot(elecps,eeg.chanlocs,'maplimits',[0,0.05],'plotchans',goods,'style','map'); colormap parula; h=colorbar; title(h,'p'); title('p < 0.05'); 

[sv,si] = sort(mdists,'descend');
subplot(4,8,5);
bar(mean(mgamma(si,:),2),'w') ; xlabel(' \leftarrow high dist (subjects) low dist \rightarrow'); ylabel('\muV gamma'); title('gamma amp. sorted by distance'); hold on; 
bar(1:12,mean(mgamma(si(1:12),:),2),'m') ;bar(24-11:24,mean(mgamma(si(end-11:end),:),2),'c') ; xlim([0,25]);
subplot(4,8,6);
bar(mean(malpha(si,:),2),'w') ; xlabel(' \leftarrow high dist (subjects) low dist \rightarrow'); ylabel('\muV alpha/beta'); title('alpha/beta amp. sorted by distance'); hold on ;
bar(1:12,mean(malpha(si(1:12),:),2),'m') ;bar(24-11:24,mean(malpha(si(end-11:end),:),2),'c') ;xlim([0,25]); 

dist_gamma = [mean(mean(gamma(si(1:12),:,elecns),2),3),mean(mean(gamma(si(end-11:end),:,elecns),3),2)];
dist_alpha = [mean(mean(alpha(si(1:12),:,elecns),2),3),mean(mean(alpha(si(end-11:end),:,elecns),3),2)];

subplot(4,8,7);
[h,p,ci,stats] = ttest2(dist_gamma(:,1),dist_gamma(:,2)); 
bar(1,mean(dist_gamma(:,1)),'m');hold on ;bar(2,mean(dist_gamma(:,2)),'c'); errorbar(mean(dist_gamma,1),std(dist_gamma,0,1)/sqrt(10),'k.') ; ylabel('\muV gamma');
title([format_t(stats.tstat),' ',format_p(p)]); set(gca,'XTick',[1,2],'XTickLabel',{'high dist','low dist'});
subplot(4,8,8); 
[h,p,ci,stats] = ttest2(dist_alpha(:,1),dist_alpha(:,2)); 
bar(1,mean(dist_alpha(:,1)),'m');hold on ;bar(2,mean(dist_alpha(:,2)),'c'); errorbar(mean(dist_alpha,1),std(dist_alpha,0,1)/sqrt(10),'k.') ; ylabel('\muV alpha/beta');
title([format_t(stats.tstat),' ',format_p(p)]) ; set(gca,'XTick',[1,2],'XTickLabel',{'high dist','low dist'});

%inot
subplot(4,8,1+8)
plot(binots',mean(mgamma,2),'kd','LineWidth',2); lsline; xlabel('I0'); ylabel('\muV gamma'); 
[rho,p] = corr(mean(mgamma,2),binots'); 
title([format_rho(rho),' ',format_p(p)]);
subplot(4,8,2+8)
plot(binots',mean(malpha,2),'kd','LineWidth',2); lsline; xlabel('I0'); ylabel('\muV alpha/beta'); 
[rho,p] = corr(mean(malpha,2),binots'); 
title([format_rho(rho),' ',format_p(p)]);

% optimal correlation threshold
% frequency vs roi size


[inotgammacorrs,inotgammaps] = corr(mean(mgamma,2),inots); 
[inotalphacorrs,inotalphaps] = corr(mean(malpha,2),inots); 

subplot(4,8,3+8);
shadedErrorBar(100:250:5000,mean(inots,1),std(inots,0,1)/sqrt(24)); xlabel('# normals'); ylabel('I0'); xlim([100,5000]) ; ylim([0.7,0.96]); title('I0 vs #normals');

subplot(4,8,4+8); 
plot(100:250:5000,inotgammacorrs,'r','LineWidth',2); hold on  ;plot(100:250:5000,inotalphacorrs,'b','LineWidth',2); hline(0,'k'); xlim([100,5000]); xlabel('# normals'); ylabel('correlation (rho)'); title('rho vs #normals');
legend({'gamma','alpha/beta'});

%{
[inotcorrs,inotps] = corr(squeeze(mean(mean(mersp(:,:,:,times>0.5 & times<2),2),4)),inots); 
subplot(4,8,3+8); 
imagesc(normls,freqs,inotcorrs); h = colorbar ; title(h,'rho'); xlabel('# voxel'); ylabel('frequency(hz)'); title('all freqs vs I0'); 
subplot(4,8,4+8) ; imagesc(normls,freqs,inotps,[0,0.1]); h=colorbar ; title(h,'p'); title('p<0.1'); xlabel('# voxel'); ylabel('frequency(hz)');
%}



[isv,isi] = sort(binots,'descend');
subplot(4,8,5+8);
bar(mean(mgamma(isi,:),2),'w') ; xlabel(' \leftarrow high I0 (subjects) low i0 \rightarrow'); ylabel('\muV gamma'); title('gamma amp. sorted by I0'); hold on; 
bar(1:12,mean(mgamma(isi(1:12),:),2),'m') ;bar(24-11:24,mean(mgamma(isi(end-11:end),:),2),'c') ; xlim([0,25]);
subplot(4,8,6+8);
bar(mean(malpha(isi,:),2),'w') ; xlabel(' \leftarrow high I0 (subjects) low i0 \rightarrow'); ylabel('\muV alpha/beta'); title('alpha/beta amp. sorted by I0'); hold on ;
bar(1:12,mean(malpha(isi(1:12),:),2),'m') ;bar(24-11:24,mean(malpha(isi(end-11:end),:),2),'c') ; xlim([0,25]);

inot_gamma = [mean(mean(gamma(isi(1:12),:,elecns),2),3),mean(mean(gamma(isi(end-11:end),:,elecns),3),2)];
inot_alpha = [mean(mean(alpha(isi(1:12),:,elecns),2),3),mean(mean(alpha(isi(end-11:end),:,elecns),3),2)];

subplot(4,8,7+8);
[h,p,ci,stats] = ttest2(inot_gamma(:,1),inot_gamma(:,2)); 
bar(1,mean(inot_gamma(:,1)),'m');hold on ;bar(2,mean(inot_gamma(:,2)),'c'); errorbar(mean(inot_gamma,1),std(inot_gamma,0,1)/sqrt(10),'k.') ; ylabel('\muV gamma');
title([format_t(stats.tstat),' ',format_p(p)]); set(gca,'XTick',[1,2],'XTickLabel',{'high I0','low I0'});
subplot(4,8,8+8); 
[h,p,ci,stats] = ttest2(inot_alpha(:,1),inot_alpha(:,2)); 
bar(1,mean(inot_alpha(:,1)),'m');hold on ;bar(2,mean(inot_alpha(:,2)),'c'); errorbar(mean(inot_alpha,1),std(inot_alpha,0,1)/sqrt(10),'k.') ; ylabel('\muV alpha/beta');
title([format_t(stats.tstat),' ',format_p(p)]) ; set(gca,'XTick',[1,2],'XTickLabel',{'high I0','low I0'});


% combining effects of i0 and distance
combs = (mat2gray(binots))/4 + (mat2gray(sqrt(mdists)))';
set(fig3,'Units','normalized');
set(fig3,'Position',[0 0 1.2 .8]);


figure,
% supplementary figure:
comp_gamma = load('allsub_gamma.mat'); comp_gamma = comp_gamma.allsub_gamma; 
comp_alpha = load('allsub_alpha.mat'); comp_alpha = comp_alpha.allsub_alpha; 
clear corrs ps ; 
for i=1:64 ; subplot(5,13,i) ; 
    plot(mdists,squeeze(mean(comp_gamma(:,i,elecns),3)),'kd','LineWidth',2) ; lsline ; 
    [corrs(i),ps(i)] = corr(mdists,squeeze(mean(comp_gamma(:,i,elecns),3)));     
end


fig4 = figure; 
subplot(4,8,1);
plot(mdists,squeeze(mean(comp_gamma(:,2,elecns),3)),'kd','LineWidth',2) ; lsline ; 
[rho,p] = corr(mdists,squeeze(mean(comp_gamma(:,2,elecns),3))); title([format_rho(rho),' ',format_p(p)]); xlabel('distance(mm)') ;ylabel('\muV gamma')
subplot(4,8,2);
plot(mdists,squeeze(mean(comp_gamma(:,5,elecns),3)),'kd','LineWidth',2) ; lsline ; 
[rho,p] = corr(mdists,squeeze(mean(comp_gamma(:,5,elecns),3)));  title([format_rho(rho),' ',format_p(p)]);xlabel('distance(mm)') ;ylabel('\muV gamma')

subplot(4,8,3); 
plot(mdists,squeeze(mean(comp_gamma(:,64,elecns),3)),'kd','LineWidth',2) ; lsline ; 
[rho,p] = corr(mdists,squeeze(mean(comp_gamma(:,64,elecns),3)));  title([format_rho(rho),' ',format_p(p)]);xlabel('distance(mm)') ;ylabel('\muV gamma')

subplot(4,8,4);
bar(corrs); ylim([-0.65,0]);xlim([0,64]); xlabel('# components'); ylabel('rho'); 

subplot(4,8,5); 
bar(ps) ; hline(0.05,'r'); ylim([0,0.3]); xlim([0,64]); xlabel('# components'); ylabel('p');

set(fig4,'Units','normalized');
set(fig4,'Position',[0 0 1.2 .8]);



pmersp = squeeze(mean(mean(mean(allstim_ersp(:,1:3,[1,5],:,times>0.5 & times<2),3),5),2)); 
for i=1:size(mersp,1)
    gamma_peak(i) = find(mersp(i,20:end)==max(mersp(i,20:end))); 
    alpha_peak(i) = find(mersp(i,4:14)==min(mersp(i,4:14))); 

end
peaksubs = [3,4,5,6,7,8,10,11,12,15,17,18,19,20,21,22,23];
for i=1:24 ; subplot(3,8,i) ; imagesc(squeeze(mean(mean(allstim_ersp(i,1:3,[1,5],:,:),3),2)),[-2,2]) ; axis xy ; colormap jet ; title(i); end

gamma_peak = (gamma_peak + 19)*2;
alpha_peak = (alpha_peak + 3)*2;

fig5=figure,
subplot(4,8,1);
plot(mdists(peaksubs),gamma_peak(peaksubs),'kd','LineWidth',2); lsline ; xlabel('distance (mm)') ; ylabel('gamma peak (Hz)'); 
[rho,p] = corr(mdists(peaksubs),gamma_peak(peaksubs)');
title([format_rho(rho),' ',format_p(p)]);
subplot(4,8,2);
plot(mdists(peaksubs),alpha_peak(peaksubs),'kd','LineWidth',2); lsline; xlabel('distance (mm)') ; ylabel('alpha/beta peak (Hz)'); 
[rho,p] = corr(mdists(peaksubs),alpha_peak(peaksubs)');
title([format_rho(rho),' ',format_p(p)]);
subplot(4,8,9);
plot(binots(peaksubs),gamma_peak(peaksubs),'kd','LineWidth',2); lsline; xlabel('I0') ; ylabel('gamma peak (Hz)'); 
[rho,p] = corr(binots(peaksubs)',gamma_peak(peaksubs)');
title([format_rho(rho),' ',format_p(p)]);
subplot(4,8,10); 
plot(binots(peaksubs),alpha_peak(peaksubs),'kd','LineWidth',2); lsline; xlabel('I0') ; ylabel('alpha/beta peak (Hz)'); 
[rho,p] = corr(binots(peaksubs)',alpha_peak(peaksubs)');
title([format_rho(rho),' ',format_p(p)]);
set(fig5,'Units','normalized');
set(fig5,'Position',[0 0 1.2 .8]);





