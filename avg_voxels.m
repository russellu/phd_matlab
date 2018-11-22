clear all ; close all 
comps = {[32,6,48],[9,18,20],[12,10,4],[24,9,16],[8,15,52],[11,21,33],[46,33,15],[25,16,21],[17,7,21]} ;

subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
subns = [1,2,3,4,5,6,7,9] ; 
for sub=1:length(subns)
    cd(['c:/shared/newbadger_mri/',subs{subns(sub)}]) ;
    rest = load_untouch_nii('singletrial_1.nii.gz') ; 
    allrests(:,:,:,:,sub) = rest.img ; 
    disp(subs{subns(sub)}) ; 


end

cd c:/shared/newbadger_mri/valerie ; 
cat = load_untouch_nii('cat.nii.gz') ; cat.img = squeeze(mean(allrests,5)) ; 
save_untouch_nii(cat,'mean_single1.nii.gz') ; 
