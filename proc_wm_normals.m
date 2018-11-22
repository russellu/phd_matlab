clear all ; close all ; 
cd C:\shared\all_white_normals\a1_good ; 
subs = dir('sub_*') ; subs(2) = []  ;

elecorder = {'FP1','FPZ','FP2','AF8','AF4','GND','AF3','AF7','F7','F5','F3','F1','FZ','F2','F4','F6','F8','FT10','FT8','FC6','FC4','FC2','REF','FC1','FC3','FC5','FT7','FT9',...
    'T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP10','TP8','CP6','CP4','CP2','CPZ','CP1','CP3','CP5','TP7','TP9','P7','P5','P3','P1','PZ','P2','P4','P6','P8',...
    'PO8','PO4','POZ','PO3','PO7','O1','OZ','O2'} ; 

goodelecs = [41,42,53,54,59,60] ; 

for s=1:length(subs) ; disp(subs(s).name) ; 
    cd(['C:\shared\all_white_normals\a1_good\',subs(s).name]) ; 
    normsleft = load_untouch_nii(['complete2whitenormals2white',subs(s).name,'_lh.mgz.nii.gz']) ; 
    normsright = load_untouch_nii(['complete2whitenormals2white',subs(s).name,'_rh.mgz.nii.gz']) ; 
    
    normsboth = normsleft.img + normsright.img ; 
    
    roi = load_untouch_nii('roi.nii.gz') ; 
    corrs = load_untouch_nii('meancorrs_t1.nii.gz') ; 
    
    roicorrs = double(corrs.img .* roi.img).*(sum(normsboth.^2,4)>0) ;
    zimg = zeros(size(roicorrs)) ; 
    inds = find(roicorrs>.1) ; 
    [sv,si] = sort(roicorrs(inds),'descend') ; 
    
    % electrode_coordinates
    coords = load('coords') ; coords = coords.coords ;
    indshigh = find(roicorrs>.3) ; 
    % for correlation threshold:
    [cx,cy,cz] = ind2sub(size(roicorrs),indshigh) ; 
    clear dipoles ; 
    for j=1:length(cx)
       dipoles(j,:) = [cx(j);cy(j);cz(j);squeeze(normsboth(cx(j),cy(j),cz(j),:))] ;       
    end  
   corrpotential = zeros(1,65) ; 
    for c=1:size(coords,2)
        for d=1:size(dipoles,1)
            dvec = coords(:,c) - dipoles(d,1:3)' ; 
            normdvec = dvec./sqrt(sum(dvec.^2)) ; 
            angle = subspace(normdvec,dipoles(d,4:6)') ; 
            corrpotential(c) = corrpotential(c) + (cos(angle)./(4*pi*sqrt(sum(dvec.^2)).^2)) ; 
        end
    end
    subcorrpots(s,:) = corrpotential ; 
    
    % number of active voxels
    numvox(s) = length(cx) ; 
    
    
    
    % for equal sized ROI
    icount =1  ;
    for i=100:800:3200
        [cx,cy,cz] = ind2sub(size(roicorrs),inds(si(1:i))) ; 
        
        % get the dipoles in the v1 mask (a dipole is a unit vector 
        % with 6 elements (position, orientation)
        clear dipoles ; 
        for j=1:length(cx)
           dipoles(j,:) = [cx(j);cy(j);cz(j);squeeze(normsboth(cx(j),cy(j),cz(j),:))] ;       
        end
        

        % calculate the potential at each electrode
        elecpotential = zeros(1,65) ; 
        for c=1:size(coords,2)
            for d=1:size(dipoles,1)
                dvec = coords(:,c) - dipoles(d,1:3)' ; 
                normdvec = dvec./sqrt(sum(dvec.^2)) ; 
                angle = subspace(normdvec,dipoles(d,4:6)') ; 
                elecpotential(c) = elecpotential(c) + (cos(angle)./(4*pi*sqrt(sum(dvec.^2)).^2)) ; 
            end
        end
        
        subelecs(s,icount,:) = elecpotential ; 
        i0s(s,icount) = 1-sqrt(sum(sum(dipoles(:,4:end),1).^2))./size(dipoles,1) ;
        icount = icount +1 ; 
        
    end
    
    elecs = load_untouch_nii('elecbrain.nii.gz') ; 
    bersp = load('bersp.mat') ; bersp = bersp.bersp ; mbersp = squeeze(mean(mean(mean(mean(bersp(1:6,1:end,:,:,25:80),1),2),3),5)) ; 
    allbersps(s,:) = mbersp ; stimbersp = squeeze(mean(mean(mean(bersp(:,1:end,:,:,25:80),2),3),5)) ; allstimbersps(s,:,:) = stimbersp ; 
    
    csizecount = 1 ; 
    for csize= 100:100:3000
        zimg = zeros(size(roicorrs)) ; 

        zimg(inds(si(1:csize))) =1 ; 
        [cx,cy,cz] = centmass3(zimg) ; 
        for i=1:size(coords,2)
            cdiffs(i,csizecount) = sqrt(sum([coords(1,i)-cx,coords(2,i)-cy,coords(3,i)-cz].^2)) ; 
        end
        csizecount = csizecount + 1 ; 
    end
    alldiffs(s,:,:) = cdiffs ; 
    %{
    nelecs = length(unique(elecs.img(find(elecs.img ~= 0)))) ; 
    clear cx cy cz
    for i=1:65
        [cx(i),cy(i),cz(i)] = centmass3(elecs.img==i) ; 
    end  
    coords = [cx;cy;cz] ; save('coords','coords') ; 
    %}
end



for i=1:size(allbersps,2)
    for j=1:size(alldiffs,2)
        for k=1:size(alldiffs,3)
            [c,p] = corr([squeeze(allbersps(:,i)),squeeze(alldiffs(:,j,k))]) ; 
            diffcorrs(i,j,k) = c(1,2) ; 
            diffps(i,j,k) = p(1,2) ; 
        end
    end
end

subplot(1,2,1) ; plot(squeeze(mean(mean(diffcorrs(:,goodelecs,:),2),3))) ; hline(0,'k') ; ylabel('correlation(r)') ; xlabel('frequency') ; set(gca,'XTick',1:5:size(powcorrs,3),'XTickLabel',(1:5:size(powcorrs,3))*2) ; 
vline([35,50],'k') ;
subplot(1,2,2) ; 
plot(squeeze(mean(mean(alldiffs(:,goodelecs,:),2),3)),squeeze(mean(allbersps(:,35:50),2)),'o') ; lsline ;  xlabel('distance(mm) from V1') ; ylabel('gamma (70-100Hz) modulation') ; 
[c,p] = corr([squeeze(mean(mean(alldiffs(:,goodelecs,:),2),3)),squeeze(mean(allbersps(:,35:50),2))]) ;
title(['rho=',num2str(c(1,2)),', p=',num2str(p(1,2))]) ; 

allbersps = squeeze(mean(allstimbersps(:,[1:6],:),2)) ;
for i=1:size(subelecs,2)
    for j=1:size(subelecs,3)
        for k=1:size(allbersps,2)
            powcorrs(i,j,k) = corr2(squeeze(subelecs(:,i,j)),squeeze(allbersps(:,k))) ; 
        end
    end
end
subplot(1,2,1),imagesc(squeeze(mean(diffcorrs(:,goodelecs,:),3))',[-.8,.8]); subplot(1,2,2) ; imagesc(squeeze(mean(powcorrs(:,goodelecs,:),1)),[-.8,.8])
    
subplot(1,2,1) ; 
plot(squeeze(mean(mean(powcorrs(:,[53,54,55,59,60,61],:),1),2))) ; hline(0,'k') ; ylabel('correlation(r)') ; xlabel('frequency') ; set(gca,'XTick',1:5:size(powcorrs,3),'XTickLabel',(1:5:size(powcorrs,3))*2) ; 
vline([35,50],'k') ; 
subplot(1,2,2) ; plot(squeeze(mean(mean(subelecs(:,:,[53,54,55,59,60,61]),2),3)),squeeze(mean(allbersps(:,35:50),2)),'o') ; lsline ;
[c,p] = corr([squeeze(mean(mean(subelecs(:,:,[53,54,55,59,60,61]),2),3)),squeeze(mean(allbersps(:,35:50),2))]) ;
title(['rho=',num2str(c(1,2)),', p=',num2str(p(1,2))]) ; 
xlabel('forward model (a.u)') ; ylabel('gamma (70-100Hz) modulation') ; 


cd c:/shared/badger_eeg/alex ; merged = pop_loadset('allstim_broad_merged.set') ;
labs = {merged.chanlocs.labels} ;
for i=1:length(labs)
    s = find(strcmpi(labs{i},elecorder)) ; 
    if ~isempty(s)
    labinds(i) = s ;  
    end
end
labinds(labinds==0) = 1 ; 
for i=1:60
  %  topos(i,:) = mean(powcorrs(:,labinds,i),1) ;
    topos(i,:) = mean(diffcorrs(i,labinds,:),3) ;

end
for i=1:60 ; subplot(6,10,i) ; topoplot(topos(i,:),merged.chanlocs,'maplimits',[-1,1]) ; title(['hz=',num2str(i*2)]) ; end ;
topoplot(topos(45,:),merged.chanlocs,'maplimits',[-.8,.8]) ;

% i0 correlations
for i=1:size(i0s,2)
    for j=1:size(allbersps,2)
        icorrs(i,j) = corr2(i0s(:,i),allbersps(:,j)) ; 
    end
end
subplot(1,2,2) ; plot(squeeze(mean(i0s(:,1:end),2)),squeeze(mean(allbersps(:,30:40),2)),'o') ; lsline ; 
[c,p] = corr([squeeze(mean(i0s(:,1:end),2)),squeeze(mean(allbersps(:,30:40),2))]) ;
title(['rho=',num2str(c(1,2)),', p=',num2str(p(1,2))]) ; ylabel('gamma modulation(60-80Hz)') ; xlabel('i0') ;
subplot(1,2,1) ; plot(squeeze(mean(icorrs(1:end,:)))) ; hline(0,'k') ; ylabel('correlation(r)') ; xlabel('frequency') ; set(gca,'XTick',1:5:size(powcorrs,3),'XTickLabel',(1:5:size(powcorrs,3))*2) ; 
vline([30,40],'k') ; 

clear corrs
for i=1:size(subcorrpots,2)
    for j=1:size(allbersps,2)
        threshcorrs(i,j) = corr2(subcorrpots(:,i),allbersps(:,j)) ; 
    end
end

% ROI size correlations
clear roicorrs
for i=1:size(allbersps,2)
   roicorrs(i) = corr2(numvox,allbersps(:,i)') ;  
end




% correlations taking distance and i0 into account:
for i=1:size(allbersps,2)
    for j=1:size(alldiffs,2)
        for k=1:size(alldiffs,3)
            bcorrs(i,j,k) = corr2(squeeze(alldiffs(:,j,k))+90*mean(i0s(:,1:end),2),allbersps(:,i)) ; 
        end
    end
end
plot(squeeze(mean(mean(alldiffs(:,goodelecs,:),3),2)).*mean(i0s(:,1:end),2),mean(allbersps(:,35:45),2),'o') ; 
title(corr2(squeeze(mean(mean(alldiffs(:,goodelecs,:),3),2)).*mean(i0s(:,1:end),2),mean(allbersps(:,35:45),2))) ; lsline 


for i=1:65;
    for j=1:30
        peakcorrs(i,j) = corr2(mean(allbersps(:,30:end),2),squeeze(alldiffs(:,i,j))) ; 
    end
end




% regression for each band
mdiffs = squeeze(mean(mean(alldiffs(:,goodelecs,:),2),3)) ; 
mi0s = squeeze(mean(i0s(:,1:2),2)) ; 
x = [ones(21,1),(mdiffs)./sqrt(sum(mdiffs.^2)),(mi0s)./sqrt(sum(mi0s.^2))] ; clear regs 
for i=1:size(allbersps,2)
   regs(i,:) = regress(squeeze(allbersps(:,i)),x) ; 
   c = corr2(allbersps(:,i),x*regs(i,:)') ; 
   subplot(6,10,i) ; 
   plot(allbersps(:,i),x*regs(i,:)','o') ; lsline ; title(c) ; 
end


plot(regs(:,2)) ; hold on ; plot(regs(:,3),'r') ; hline(0,'k')










