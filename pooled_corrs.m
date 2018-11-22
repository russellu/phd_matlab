clear all ;close all ; 
subs = {'a_alex','a_charest','a_esteban','a_fabio','a_gab','a_gabriella','a_genevieve','a_gina','a_guillaume','a_jeremie','a_julie','a_katrine','a_lisa','a_marc',...
    'a_marie','a_mathieu','a_maxime','a_mingham','a_patricia','a_po','a_russell','a_sunachakan','a_tah','a_vincent',...
    'b_alex','b_dina','b_genevieve','b_jeremie','b_karl','b_russell','b_sukhman','b_tegan','b_valerie'} ; 

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
    [sv,si] = sort(corrs.img(:),'descend'); 
    zcorrs(si(1:10000)) = 1; 
    
    % get distance from ROI
    [cx,cy,cz] = centmass3(zcorrs); 
    coord_diffs = [coords(:,1) - cx, coords(:,2) - cy, coords(:,3) - cz];
    dists = sqrt(sum(coord_diffs.^2,2)); 
    subdists(sb,:) = dists; 

end



cd E:\nimg_pool\saved
a_elabs = load('a_elabs'); a_elabs = a_elabs.a_elabs; 
badger_elabs = load('badger_elabs'); badger_elabs = badger_elabs.badger_elabs; 
elecorder = load('elecorder'); elecorder = elecorder.elecorder; 
all_postchans = {'O1','OZ','OZ','PO3','POZ','PO4','P8','P6','P4','P2','PZ','P1','P3','P5','P7'}; % 
for i=1:length(all_postchans)
   a_chaninds(i) = find(strcmpi(all_postchans{i},a_elabs)); 
   b_chaninds(i) = find(strcmpi(all_postchans{i},badger_elabs)); 
   elecorder_chaninds(i) = find(strcmpi(all_postchans{i},elecorder)); 
end

allsub_gamma = load('allsub_gamma'); allsub_gamma = allsub_gamma.allsub_gamma; 
allsub_gamma_9 = load('allsub_gamma_9'); allsub_gamma_9 = allsub_gamma_9.allsub_gamma_9; 

bothsub_gamma = zeros(33,64,15); bothsub_gamma(1:24,:,:) = allsub_gamma(:,:,a_chaninds); bothsub_gamma(25:end,:,:) = allsub_gamma_9(:,:,b_chaninds); 
bothsub_dists = zeros(1,33) ; bothsub_dists(1:24) = mean(subdists(1:24,elecorder_chaninds),2); bothsub_dists(25:end) = mean(subdists(25:end,elecorder_chaninds),2); 
subinds = 1:24; 
for i=1:64 ; [rhos(i),ps(i)] = corr(bothsub_dists(subinds)',squeeze(mean(bothsub_gamma(subinds,i,:),3))); end

subplot(1,5,1);
bar(rhos); xlabel('# components'); ylabel('spearman rho') ; xlim([1,64])
subplot(1,5,2);
bar(ps) ; hline(0.05,'r'); xlabel('# components'); ylabel('p-value'); xlim([1,64]);
subplot(1,5,3);
plot(bothsub_dists(subinds),squeeze(mean(bothsub_gamma(subinds,2,:),3)),'kd','LineWidth',2); title(['rho=',num2str(rhos(2)),' p=',num2str(ps(2))]); xlabel('distance(mm)'); ylabel('gamma amplitude (micro-volt)'); lsline
subplot(1,5,4);
plot(bothsub_dists(subinds),squeeze(mean(bothsub_gamma(subinds,5,:),3)),'kd','LineWidth',2); title(['rho=',num2str(rhos(5)),' p=',num2str(ps(5))]); xlabel('distance(mm)'); ylabel('gamma amplitude (micro-volt)'); lsline
subplot(1,5,5);
plot(bothsub_dists(subinds),squeeze(mean(bothsub_gamma(subinds,64,:),3)),'kd','LineWidth',2); title(['rho=',num2str(rhos(64)),' p=',num2str(ps(64))]); xlabel('distance(mm)'); ylabel('gamma amplitude (micro-volt)'); lsline




