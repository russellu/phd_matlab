clear all ;close all ; 
subs = {'a_alex','a_charest','a_esteban','a_fabio','a_gab','a_gabriella','a_genevieve','a_gina','a_guillaume','a_jeremie','a_julie','a_katrine','a_lisa','a_marc',...
    'a_marie','a_mathieu','a_maxime','a_mingham','a_patricia','a_po','a_russell','a_sunachakan','a_tah','a_vincent'} ; 

eeg = pop_loadset('e:/nimg_pool/a_alex/cleanfilt.set') ; 


for sb=1:length(subs)
    
    % DISTANCE
    cd(['e:/nimg_pool/',subs{sb}]) ; 
    alphabeta = load_untouch_nii('alphabeta_fs.nii.gz');
    gamma = load_untouch_nii('gamma_fs.nii.gz'); 
    corrs = load_untouch_nii('cleancorrs_fs.nii.gz'); 
    corrs.img(:,:,100:end) = 0; 
    zcorrs = zeros(size(corrs.img)); 
    corrs.img(isnan(corrs.img)) = 0; 
    [sv,si] = sort(corrs.img(:),'descend'); 
    
    [gsv,gsi] = sort(gamma.img(:),'descend'); 
    [asv,asi] = sort(alphabeta.img(:),'ascend'); 
    
    zgamma = zeros(size(gamma.img));
    zalphabeta = zeros(size(alphabeta.img)); 
    
    nvox = 15000; 
    
    zgamma(gsi(1:nvox)) = 1;
    zalphabeta(asi(1:nvox)) = 1;
    zcorrs(si(1:nvox)) = 1; 
   %figure,subplot(2,2,1) ; imagesc(squeeze(mean(zcorrs,3))); 
   % subplot(2,2,2) ; imagesc(squeeze(mean(zgamma,3))); 
   % subplot(2,2,3) ; imagesc(squeeze(mean(zalphabeta,3))); 
    
    src_gamma(sb) = mean(gamma.img(si(1:nvox))); 
    src_alphabeta(sb) = mean(alphabeta.img(si(1:nvox))); 
    
    olapgamma(sb) = sum(sum(sum(zgamma.*zcorrs)))/nvox;
    olapalpha(sb) = sum(sum(sum(zalphabeta.*zcorrs)))/nvox; 

    coords = load('coords.mat'); coords = coords.coords; 
    size(coords); 
    subcoords(sb,:,:) = coords; 
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
    
    allgammas(sb,:,:,:) = gamma.img; 
    
end

for i=1:size(allgammas,1)
    gammai = squeeze(allgammas(i,:,:,:)); 
    gammai(isnan(gammai)) = 0; 
    gammai(:,:,100:end) = 0; 

    [sv,si] = sort(gammai(:),'descend'); 
    zgammai = zeros(size(gammai));
    zgammai(si(1:15000)) = 1; 
    [cx,cy,cz] = centmass3(zgammai); 
    
    coord_diffs = [coords(:,1) - cx, coords(:,2) - cy, coords(:,3) - cz];
    dists = sqrt(sum(coord_diffs.^2,2)); 
    subsrcdists(i,:) = dists; 
    
    roisrcgamma(i) =  mean(gammai(si(1:15000)))/1000000000; 

end

cd E:\nimg_pool\saved
allstim_ersp = load('allstim_ersp'); allstim_ersp = allstim_ersp.allstim_ersp; 
mersp = squeeze(mean(allstim_ersp(:,1:6,:,:,:),2)); 
times = load('times'); times = times.times ; freqs = load('freqs') ; freqs = freqs.freqs; 

plot(squeeze(mean(mean(mean(mersp(:,[1,5],freqs>40 & freqs<90,times>0 & times<2),2),3),4)),src_gamma,'kd','LineWidth',2) ; lsline
[rho,p] = corr(squeeze(mean(mean(mean(mersp(:,[1,5],freqs>40 & freqs<90,times>0 & times<2),2),3),4)),src_gamma'); 


cd E:\nimg_pool\saved
times = load('times'); times = times.times ; freqs = load('freqs') ; freqs = freqs.freqs; 
a_elabs = load('a_elabs'); a_elabs = a_elabs.a_elabs; 
elecorder = load('elecorder'); elecorder = elecorder.elecorder; 
all_postchans = {'O1','OZ','OZ','PO3','POZ','PO4','P8','P6','P4','P2','PZ','P1','P3','P5','P7'}; %
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

mdists = squeeze(mean(subdists(:,elecorder_chaninds),2)); 

msrcdists = squeeze(mean(subsrcdists(:,elecorder_chaninds),2)); 



allsub_allstim_alpha = load('allsub_allstim_alpha'); alpha = allsub_allstim_alpha.allsub_allstim_alpha; 
allsub_allstim_gamma = load('allsub_allstim_gamma'); gamma = allsub_allstim_gamma.allsub_allstim_gamma; 
allstim_ersp = load('allstim_ersp'); allstim_ersp = allstim_ersp.allstim_ersp; 
mersp = squeeze(mean(allstim_ersp(:,1:6,:,:,:),2)); 

mgamma = squeeze(mean(gamma(:,:,a_chaninds),3)); 
malpha = squeeze(mean(alpha(:,:,a_chaninds),3)); 


src_gamma = src_gamma/1000000000;
src_alphabeta = src_alphabeta/1000000000; 

fig = figure;
subplot(4,8,1);
plot(squeeze(mean(mgamma(:,:),2)),src_gamma,'kd','LineWidth',2) ; lsline ;  xlabel('scalp gamma amp. (\muV)'); ylabel('source gamma amp. (in BOLD ROI) (A.U)')
[rho,p] = corr(src_gamma',squeeze(mean(mgamma(:,:),2)));
title([format_rho(rho),' ',format_p(p)]);
subplot(4,8,3);
plot(mdists,src_gamma,'kd','LineWidth',2); lsline; xlabel('BOLD ROI distance (mm)'); ylabel('source gamma amp. (in BOLD ROI) (A.U)'); 
[rho,p] = corr(mdists,src_gamma');
title([format_rho(rho),' ',format_p(p)]);

subplot(4,8,2); 
plot(squeeze(mean(malpha(:,:),2)),src_alphabeta,'kd','LineWidth',2) ; lsline ;xlabel('scalp alpha/beta amp. (\muV)'); ylabel('source alpha/beta amp. (in BOLD ROI) (A.U)')
[rho,p] = corr(squeeze(mean(malpha(:,:),2)),src_alphabeta');
title([format_rho(rho),' ',format_p(p)]);
subplot(4,8,4);
plot(mdists,src_alphabeta,'kd','LineWidth',2); lsline; xlabel('BOLD ROI distance (mm)'); ylabel('source alpha/beta amp. (in BOLD ROI) (A.U)'); 
[rho,p] = corr(mdists,src_alphabeta');
title([format_rho(rho),' ',format_p(p)]);

subplot(4,8,5); 
olaps = [olapgamma;olapalpha]*100;
barwitherr(squeeze(std(olaps,0,2))/sqrt(24),squeeze(mean(olaps,2))); 
[h,p,ci,stats] = ttest(olaps(1,:),olaps(2,:)); 
title([format_t(stats.tstat),' ',format_p(p)]); set(gca,'XTickLabel',{'gamma','alpha/beta'}); ylabel('overlap with BOLD ROI (%)'); xlabel('ROI'); 

subplot(4,8,6); 
plot(msrcdists,roisrcgamma','kd','LineWidth',2); lsline;  xlabel('source ROI distance (mm)'); ylabel('source gamma amp. (in source ROI) (A.U)'); 
[rho,p] = corr(msrcdists,src_gamma'); 
title([format_rho(rho),' ',format_p(p)]);

subplot(4,8,7); 
plot(msrcdists,mdists,'kd','LineWidth',2); lsline;  xlabel('source ROI distance (mm)'); ylabel('BOLD ROI distance (mm)'); 
[rho,p] = corr(msrcdists,mdists); 
title([format_rho(rho),' ',format_p(p)]);



[sv,si] = sort(msrcdists,'descend'); 
sorted_roisrcgamma = [roisrcgamma(si(1:12));roisrcgamma(si(13:end))]; 
sorted_msrcdists = [msrcdists(si(1:12))';msrcdists(si(1:12))'];

subplot(4,8,8); 
bar(mean(sorted_roisrcgamma,2)) ; hold on ; errorbar(squeeze(mean(sorted_roisrcgamma,2)),std(sorted_roisrcgamma,0,2)/sqrt(12),'k.'); 
[h,p,ci,stats] = ttest2(sorted_roisrcgamma(1,:),sorted_roisrcgamma(2,:)); set(gca,'XTickLabel',{'high dist','low dist'}); ylabel('source gamma amp. (in source ROI) (A.U)'); xlabel('source to electrode dist');
title([format_t(stats.tstat),' ',format_p(p)]);

set(fig,'Units','normalized');
set(fig,'Position',[0 0 1.2 .8]);