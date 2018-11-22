clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 
threshs = 0.1:0.05:0.8; 
for sb=1:length(subs); 
    cd(['c:/shared/allfmris/sub_',subs{sb}]) ;  
    corrs = load_untouch_nii('cleancorrs_fs.nii.gz'); 
    corrcopy = corrs; 
    [sv,si] = sort(corrs.img(:),'descend') ; 
    zcorrs = zeros(size(corrs.img)) ; zcorrs(si(1:8000)) = 1 ; 
    [sx,sy,sz] = ind2sub(size(corrs.img),si(1:8000)) ; 
    [cx,cy,cz] = centmass3(zcorrs) ; 
    
    zimg = zeros(size(corrs.img)); 
    zimg(cx-3:cx+3,cy-3:cy+3,cz-3:cz+3) = 1; 
    
    %[gx,gy,gz] = ndgrid(-size(zcorrs,1)/2:size(zcorrs,1)/2-1,-size(zcorrs,2)/2:size(zcorrs,2)/2-1,-size(zcorrs,3)/2:size(zcorrs,3)/2-1);
    %gx = gx+cx; gy = gy+cy; gz = gz+cz; 
    %sphere = sqrt(gx.^2 + gy.^2 + gz.^2) > 15; 
    
    corrcopy.img = zimg; save_untouch_nii(corrcopy,'centsqr.nii.gz'); 
   
    mkdir corrdir
    cd corrdir
    for thresh=1:length(threshs)
       binimg = corrs.img>threshs(thresh); 
       corrcopy.img = binimg ; save_untouch_nii(corrcopy,['thresh_',num2str(thresh),'.nii.gz']); 
    end

end