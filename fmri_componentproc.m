clear all ; close all ; 
subs = {'sub_alex','sub_charest','sub_esteban','sub_fabio','sub_gab','sub_gabriella','sub_genevieve','sub_gina','sub_jeremie','sub_julie','sub_katrine','sub_lisa'...
        ,'sub_marc','sub_marie','sub_mathieu','sub_maxime','sub_mingham','sub_patricia','sub_po','sub_russell','sub_sunachakan','sub_vincent'} ;
cd c:/shared/fmricomps ; ls
subics = load('subics') ; subics = subics.subics ; 
submaps = load('submaps') ; submaps = submaps.submaps ; 
subparams = load('subparams') ; subparams = subparams.subparams ; 
subtrigs = load('subtrigs') ; subtrigs = subtrigs.subtrigs ; 
for i=1:size(subics,1)
    ic_i = subics(i,:) ;
    scount = 1 ; 
    for sc=1:245:2205 ; scans(scount,:) = ic_i(sc:sc+244) ; scount = scount + 1 ; end
    params_i = squeeze(subparams(i,:,:)) ; 
    trigs_i = round(squeeze(subtrigs(i,:,:))/2) ; 
    stimcounts = ones(1,6) ; 
    ecount = 1 ; 
    for sc=1:size(params_i,1) ;
        for stim=1:6  
            inds = find(params_i(sc,:)==stim) ; 
            for pres=1:length(inds)
                epoch = scans(sc,(trigs_i(sc,inds(pres))-1):(trigs_i(sc,inds(pres))+8)) ; 
                stims(stim,stimcounts(stim),:) = epoch ; 
                stimcounts(stim) = stimcounts(stim) + 1 ;
            end
        end
    end
   % figure,
   % for st=1:6 ; subplot(2,3,st) ; imagesc(squeeze(stims(st,:,:)),[-5,5]) ; title(subs{i}) ; end
    allstims(i,:,:,:) = stims ; 
end
basestims = zeros(size(allstims)) ; 
for i=1:22 ; 
    for j=1:6
        for k=1:45
            basestims(i,j,k,:) = (allstims(i,j,k,:)-repmat(allstims(i,j,k,1),[1,1,1,10])) ; 
        end
    end
end
for i=1:22 ; 
    subplot(5,5,i) ; errorbar(squeeze(mean(basestims(i,[1,3,2],:,:),3))',squeeze(std(basestims(i,[1,3,2],:,:),0,3))'./sqrt(45),'LineWidth',2) ; title(subs{i}) ; 
    ylim([-.5,2.5]) ;xlim([1,10]) ; hline(0,'k') ; 
end





