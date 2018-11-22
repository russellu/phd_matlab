clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 

for sb=1:length(subs); disp(sb) ; 
    cd(['c:/shared/allfmris/sub_',subs{sb}]) ;  
    corrs = load_untouch_nii('cleancorrs_fs.nii.gz') ; 
    bincorrs = double(corrs.img>.3) ; 
    cd(['c:/shared/gamma_t1s/sub_',subs{sb},'/mri']) ;  
    leftnorms = load_untouch_nii('surf.nii.gz') ; 
    rightnorms = load_untouch_nii('surfr.nii.gz') ; 
    bothnorms = leftnorms.img + rightnorms.img ; 
    rlat = load_untouch_nii('r_latvis.nii.gz') ; 
    llat = load_untouch_nii('l_latvis.nii.gz') ; 
    
    rcalc = load_untouch_nii('r_calc.nii.gz') ; 
    lcalc = load_untouch_nii('l_calc.nii.gz') ; 
    
    latvis = double(imdilate(rlat.img+llat.img,strel(ones(3,3,3)))>0) ; 
    calc = double(imdilate(rcalc.img+lcalc.img,strel(ones(3,3,3)))>0) ; 

    latvisperc(sb) = sum(bincorrs(:).*latvis(:))./sum(latvis(:)) ; 
    calcperc(sb) = sum(bincorrs(:).*calc(:))./sum(calc(:)) ; 

    binnorms = sum(bothnorms,4)~=0 ; 
    bininds = find(binnorms.*latvis) ; 
    [cx,cy,cz] = ind2sub(size(binnorms),bininds) ; clear allnorms 
    for i=1:length(cx) ; allnorms(i,:) = squeeze(bothnorms(cx(i),cy(i),cz(i),:)) ; end
    inots(sb) = 1 - sum(abs(sum(allnorms,1)))./length(cx) ;
end
cd c:/shared/savedata_struct ; 
large_latvis_inots = inots ; save('large_latvis_inots','large_latvis_inots') ; 





