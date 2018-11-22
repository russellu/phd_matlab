clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 

for sb=1%:length(subs); disp(sb) 
    
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
        
        zinds = zeros(size(zcorrs)); 
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
            zinds(inds) = 1 ;

        end
        brain.img = zinds ; 
        save_untouch_nii(brain,'zinds.nii.gz'); 
end

