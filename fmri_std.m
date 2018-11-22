clear all ; close all ; 
subs = {'sub_alex','sub_charest','sub_esteban','sub_fabio','sub_gab','sub_gabriella','sub_genevieve','sub_gina','sub_jeremie','sub_julie','sub_katrine','sub_lisa'...
        ,'sub_marc','sub_marie','sub_mathieu','sub_maxime','sub_mingham','sub_patricia','sub_po','sub_russell','sub_sunachakan','sub_vincent'} ;
    
for sub=1:length(subs)
    cd(['c:/shared/allfmris/',subs{sub}]) ; ls
    cat = load_untouch_nii('cat.nii.gz') ;
    f1 = load_untouch_nii('f1.nii.gz') ; 
    catim = cat.img ; 
    f1.img = (squeeze(mean(catim,4))./squeeze(std(catim,0,4))) ;
    save_untouch_nii(f1,'coeff_var.nii.gz'); 
    f1.img = (squeeze(std(catim,0,4))) ;
    save_untouch_nii(f1,'std.nii.gz'); 

    
    
end