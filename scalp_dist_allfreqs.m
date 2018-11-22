clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 

for sb=1:length(subs); disp(sb) 
    cd(['c:/shared/freesurfer_segs/sub_',subs{sb},'/mri']) ;  
    ctx = load_untouch_nii('cortex.nii.gz');    
    %[cort_x,cort_y,cort_z] = ind2sub(size(ctx.img),find(ctx.img>0));    
    ctximg = double(imdilate(ctx.img,strel(ones(5,5,5)))); 

    cd(['c:/shared/allfmris/sub_',subs{sb}]) ;  
    brain = load_untouch_nii('fs_brain.nii.gz'); 
    tcorrs = load_untouch_nii('cleancorrs_fs.nii.gz'); 
    tcorrs.img(:,:,100:end) = 0; 
    %tcorrs.img = imfilter(double(tcorrs.img),fspecial('gaussian',11,9)); 
    [sv,si] = sort(tcorrs.img(:),'descend'); 
    zcorrs = zeros(size(tcorrs.img));
    zcorrs(si(1:10000)) = 1; 
    zcorrs = double(zcorrs).*tcorrs.img; 
    
    %[cort_x,cort_y,cort_z] = ind2sub(size(zcorrs),find(zcorrs>0)); 
    
    [cort_x,cort_y,cort_z] =   ind2sub(size(zcorrs),find(zcorrs>0));
    [cort_x2,cort_y2,cort_z2] =  centmass3((zcorrs)); 

       
    locs = load('locs.mat'); 
    locs = locs.locs; 
    
    for i=1:size(locs,2)
        dists_i = sqrt((cort_x-locs(1,i)).^2 + (cort_y-locs(2,i)).^2 + (cort_z-locs(3,i)).^2); 
        mindists(sb,i) = min(dists_i); 
        comdists(sb,i) = sqrt((cort_x2-locs(1,i)).^2 + (cort_y2-locs(2,i)).^2 + (cort_z2-locs(3,i)).^2); 
    end
    cd(['c:/shared/allres/',subs{sb}]) ;  
    mersp = load('amersp') ; mersp = mersp.amersp ; allmersp(sb,:,:,:) = mersp ;     
        
     cd(['E:\clean_allres\',subs{sb}]) ; ls 
     freqtopos = load('freqtopos.mat'); 
     all_freqtopos(sb,:,:) = freqtopos.freqtopos; 
    
end
cd(['E:\clean_allres\vincent']) ; ls 
alltopos = load('alltopos'); alltopos = alltopos.alltopos; 

cd C:\shared\allres\alex
eeg = pop_loadset('resamp_vis09.set'); 
elabs = {eeg.chanlocs.labels}; 
cd c:/shared/ ; 
elecorder = load('elecorder.mat'); elecorder = elecorder.elecorder; 
for i=1:length(elecorder)
    eleci = elecorder{i};
    indi = find(strcmpi(eleci,elabs)); 
    if ~isempty(indi)
        sbothdists(:,indi) = comdists(:,i); 
        
    end
end





