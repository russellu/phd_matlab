cd E:\ProbAtlas_v4\subj_vol_all

leftatlas = load_untouch_nii('maxprob_vol_lh.nii.gz'); 
rightatlas = load_untouch_nii('maxprob_vol_rh.nii.gz'); 
roinames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','VO1','VO2','PHC1','PHC2','MST','hMT','LO2','LO1','V3b','V3a','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5','SPL1','FEF'};
roi_inds = {[1,2],[3,4],[5,6],[7],[8,9],[10,11],[12],[13],[14,15]};
rnames = {'V1','V2','V3','hV4','VO','PH','MST','hMT','LO'};
newrois = zeros(size(leftatlas.img)); 
resnew = newrois(:); 
clear atlasinds
for i=1:length(roi_inds)
    if length(roi_inds{i}) == 2
        atlasinds{i} = [find(leftatlas.img + rightatlas.img == roi_inds{i}(1));find(leftatlas.img + rightatlas.img == roi_inds{i}(2))]; 
        resnew(atlasinds{i}) = i; 
    else
        atlasinds{i} = find(leftatlas.img + rightatlas.img == roi_inds{i}(1)); 
        resnew(atlasinds{i}) = i; 
    end
end

newrois = reshape(resnew,size(newrois)); 
leftatlas.img = newrois ; save_untouch_nii(leftatlas,'newrois.nii.gz'); 