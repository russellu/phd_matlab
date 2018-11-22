clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 

for sb=1:length(subs); disp(sb) 
    
        cd(['c:/shared/freesurfer_segs/sub_',subs{sb},'/mri']) ;  
        ctx = load_untouch_nii('cortex.nii.gz');
        ctximg = double(imdilate(ctx.img,strel(ones(5,5,5)))); 
        
        cd(['c:/shared/allfmris/sub_',subs{sb}]) ;  
        brain = load_untouch_nii('fs_brain.nii.gz'); 
        tcorrs = load_untouch_nii('cleancorrs_fs.nii.gz'); 
        tcorrs.img(:,:,100:end) = 0; 
        [sv,si] = sort(tcorrs.img(:),'descend'); 
        zcorrs = zeros(size(tcorrs.img));
        zcorrs(si(1:10000)) = 1; 
        [cx,cy,cz] = centmass3(brain.img); 
        
        subplot(4,6,sb); imagesc(squeeze(mean(zcorrs,2))); 
        locs = load('locs.mat'); 
        locs = locs.locs; 
        
        % create lines between every electrode and the ROI
        for i=1:length(locs); disp(['elec = ',num2str(i)]); 
            vecdiff = [cx - locs(1,i), cy - locs(2,i), cz - locs(3,i)]; 
            dist = round(norm(vecdiff)); 
            normdiff = vecdiff/norm(vecdiff); 
            xstep = vecdiff(1)*normdiff(1); 
            vecinds = zeros(3,dist); 
            for d=1:dist
                vecinds(1,d) = locs(1,i) + normdiff(1)*d; 
                vecinds(2,d) = locs(2,i) + normdiff(2)*d; 
                vecinds(3,d) = locs(3,i) + normdiff(3)*d; 
            end
            
            inds = sub2ind(size(zcorrs),round(vecinds(1,:)),round(vecinds(2,:)),round(vecinds(3,:))); 
            zinds = zeros(size(zcorrs)); 
            zinds(inds) = 1 ;
            maskctx =  imdilate(zinds.*ctximg,strel(ones(3,3,3))); 
            ctxinds = find(maskctx==1); 
            [cort_x,cort_y,cort_z] = ind2sub(size(maskctx),ctxinds); 
            
            elec_cort_dist = [locs(1,i) - cort_x , locs(2,i) - cort_y , locs(3,i) - cort_z]; 
            norm_elec_cort_dist = sqrt(sum(elec_cort_dist.^2,2))' ;
            mindists(sb,i) = min(norm_elec_cort_dist); 
            minbraindists(sb,i) = dist - mindists(sb,i); 
            minbothdists(sb,i) = dist; 
            
        end
        
         cd(['c:/shared/allres/',subs{sb}]) ;  
        mersp = load('amersp') ; mersp = mersp.amersp ; allmersp(sb,:,:,:) = mersp ; 
end


fmersp = squeeze(mean(mean(allmersp(:,:,:,50:180),4),2)); 
for i=1:65 
    for j=1:60
            corrs(i,j) = corr2(minbothdists(:,i),fmersp(:,j)); 
    end
end

for i=1:65
    braincorrs(i) = corr2(squeeze(mean(fmersp(:,20:50),2)),minbraindists(:,i)); 
    scalpcorrs(i) = corr2(squeeze(mean(fmersp(:,20:50),2)),mindists(:,i)); 
    bothcorrs(i) = corr2(squeeze(mean(fmersp(:,20:50),2)),minbothdists(:,i)); 

end


cd c:/shared/ ; 
elecorder = load('elecorder.mat'); elecorder = elecorder.elecorder; 
%elabs = load('elabs.mat'); elabs = elabs.elabs; 

cd C:\shared\allres\alex
eeg = pop_loadset('resamp_vis09.set'); 
elabs = {eeg.chanlocs.labels}; 

for i=1:length(elecorder)
    eleci = elecorder{i};
    indi = find(strcmpi(eleci,elabs)); 
    if ~isempty(indi)
        sbothcorrs(indi) = bothcorrs(i); 
    end
end
topoplot(sbothcorrs,eeg.chanlocs)


