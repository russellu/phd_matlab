clear all ; close all ; 
subs= {'alex','dina','genevieve','jeremie','karl','russell','tegan','valerie'} ;

for s=1:length(subs) ; 
    clear ts tsinds corrinds maxcorrs
    cd(['c:/shared/badger_mri/',subs{s},'/nii']) ; ls 
    allangles = load_untouch_nii('allangles.nii.gz') ; 
    maxangles = max(allangles.img,[],4) ; 
    [cx,cy,cz]  = ind2sub(size(maxangles),find(maxangles~=0)) ; 
    maxinds = zeros(size(maxangles)) ; 
    for i=1:length(cx) 
        maxinds(cx(i),cy(i),cz(i)) = find(squeeze(allangles.img(cx(i),cy(i),cz(i),:)) == max(squeeze(allangles.img(cx(i),cy(i),cz(i),:)))) ; 
    end
    fref = load_untouch_nii('fref.nii.gz') ; 
    fref.img = maxinds; save_untouch_nii(fref,'maxinds.nii.gz') ; 
    
end