clear all ;close all ; 
subs = {'b_alex','b_dina','b_genevieve','b_jeremie','b_karl','b_russell','b_sukhman','b_tegan','b_valerie'} ; 

eeg = pop_loadset('e:/nimg_pool/b_alex/cleanfilt.set') ; 

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

    
    %brain = load_untouch_nii('fs_t1.nii.gz'); 
    %zmask  = double(zcorrs> 0 | corrsphere.img>0); 
    %maskcorrs = zmask; maskcorrs(sphere>0) = -1; 
    %subplottight(2,12,sb);plotoverlayIntensity2D(squeeze(brain.img(:,cy,:)),squeeze(zmask(:,cy,:)),squeeze(maskcorrs(:,cy,:)),90); 

    leftnorms = load_untouch_nii('left_norms.nii.gz') ; 
    rightnorms = load_untouch_nii('right_norms.nii.gz') ; 
    bothnorms = leftnorms.img + rightnorms.img ; 
    
    norm_mask = sum(bothnorms,4)>0; 
    norm_corrs = corrs.img.*norm_mask;
    norm_corrs(isnan(norm_corrs)) = 0; 
    [sv,si] = sort(norm_corrs(:),'descend'); 
    
    normls = 100:500:8000 ; 
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


cd e:/nimg_pool/saved ; ls 
allersp = load('allersp_9'); allersp = allersp.allersp_9; 

times = load('times_9') ; times = times.times_9; 
freqs = load('freqs_9') ; freqs = freqs.freqs_9; 

cd E:\nimg_pool\saved
%badger_elabs = load('badger_elabs'); badger_elabs = badger_elabs.badger_elabs; 
badger_elabs = {eeg.chanlocs.labels};
elecorder = load('elecorder'); elecorder = elecorder.elecorder; 
mgamma = load('mgamma_9'); mgamma = mgamma.mgamma_9;
malpha = load('malpha_9'); malpha = malpha.malpha_9; 


egamma = squeeze(mean(mean(mgamma(:,:,elecns),2),3)); 

% get the indices of the labeled electrodes that match the eeg structure
for i=1:length(badger_elabs)
   indi = find(strcmpi(badger_elabs{i},elecorder));
   if ~isempty(indi)
      inds(i) = indi;         
   else
       inds(i) = 1; 
    end
end

mg = squeeze(mean(mgamma,2));    
ma = squeeze(mean(malpha,2)); 

[rhos,ps] = corr(squeeze(mean(mgamma,2)),subdists(:,inds)); 
clear rhos ps; 
for i=1:length(inds)

   [rhos(i),ps(i)] = corr(squeeze(mg(:,i)),subdists(:,inds(i)));     
end

[rho,p] = corr(subdists(:,inds(19)),mg(:,19)); 
plot(subdists(:,inds(19)),mg(:,19),'kd','LineWidth',2) ; lsline
title([format_rho(rho),' ',format_p(p)]);



fig= figure;
subplot(4,8,1); 
imagesc(times,freqs,squeeze(mean(allersp(:,1,:,:))),[-3,3]) ;axis xy ; colormap jet; xlabel('time(s)') ;ylabel('frequency(hz)'); title('0%random'); 
subplot(4,8,2); 
imagesc(times,freqs,squeeze(mean(allersp(:,2,:,:))),[-3,3]) ;axis xy ; colormap jet;  title('10%random'); 
subplot(4,8,3); 
imagesc(times,freqs,squeeze(mean(allersp(:,3,:,:))),[-3,3]) ;axis xy ; colormap jet;  title('100%random'); 
subplot(4,8,4); 
[rho,p] = corr(mean(subdists(:,inds(elecns)),2),mean(mg(:,elecns),2)); 
plot(mean(subdists(:,inds(elecns)),2),mean(mg(:,elecns),2),'kd','LineWidth',2) ; lsline ; xlabel('distance(mm)'); ylabel('\muV gamma'); 
title([format_rho(rho),' ',format_p(p)]);
subplot(4,8,5); 
[rho,p] = corr(mean(subdists(:,inds(elecns)),2),mean(ma(:,elecns),2)); 
plot(mean(subdists(:,inds(elecns)),2),mean(ma(:,elecns),2),'kd','LineWidth',2) ; lsline; xlabel('distance(mm)'); ylabel('\muV alpha/beta'); 
title([format_rho(rho),' ',format_p(p)]);
subplot(4,8,6); 
topoplot(rhos,eeg.chanlocs,'maplimits',[-.8,.8],'style','map','plotchans',goodchans); title('gamma vs distance'); colormap parula ; h = colorbar ; title(h,'rho'); 

subplot(4,8,7) ; imagesc([-3,3]) ; h = colorbar ; title(h,'dB'); 
set(fig,'Units','normalized');
set(fig,'Position',[0 0 1.2 .8]);








