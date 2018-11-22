clear all ; close all ; 

subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 
scans = {'mc_retino_allstims_01','mc_retino_allstims_02','mc_retino_gamma_01','mc_retino_gamma_02','mc_retino_movie','mc_retino_rest'}; 
eegscans = {'retino*allstim*01*set','retino*allstim*02*set','retino*gamma*01*set','retino*gamma*02*set','retino*movie*set','retino*rest*set'}; 
clean_fmri_names = {'bp_clean_retino_allstims_01','bp_clean_retino_allstims_02','bp_clean_retino_gamma_01','bp_clean_retino_gamma_02','bp_clean_retino_movie','bp_clean_retino_rest'};

mask = load_untouch_nii('c:/shared/epireg/meanmask.nii.gz'); 
maskinds = find(mask.img==1); 

xcorrs = zeros(8,6,35994,41); 
for sb=1%:length(subs)
    for sc=1%:length(scans)
        cd(['E:\rawbadger\badger_mri\',subs{sb},'\nii\',scans{sc},'.nii.mat']);
        mats = dir('MAT*');
        ls 
        
        clear allmats; 
        for mat=1:length(mats)
            m1 = load(mats(mat).name);
            allmats(:,:,mat) = m1;    
        end

        clear newmat; 
        newmat(1,:) = squeeze(allmats(1,2,:)); 
        newmat(2,:) = squeeze(allmats(1,3,:)); 
        newmat(3,:) = squeeze(allmats(2,1,:)); 
        newmat(4,:) = squeeze(allmats(2,3,:)); 
        newmat(5,:) = squeeze(allmats(3,1,:)); 
        newmat(6,:) = squeeze(allmats(3,2,:)); 

        zmat = zscore(newmat,[],2); 

        summat = sum(abs(diff(zmat,1,2))); 
        zsummat = zscore(summat); 
        zsummat = smooth(zsummat) ;
        
        cd(['e:/fmris/badger_',subs{sb},'/atlas_fmri']);
        fmri_nii = load_untouch_nii([clean_fmri_names{sc},'.nii.gz']);
        fmri_nii.img = fmri_nii.img(:,:,:,2:end); 
        
        res_fmri = reshape(fmri_nii.img,[numel(fmri_nii.img(:,:,:,1)),size(fmri_nii.img,4)]); 
        maskvox = res_fmri(maskinds,:); 
        
        for i=1:size(maskvox,1)
           xcorrs(sb,sc,i,:) = xcorr(maskvox(i,20:end-20),zsummat(20:end-20),20,'coeff');             
        end
        
        disp(subs{sb}); 
    end
end

%for i=1:6 ; subplot(2,3,i) ; imagesc(squeeze(mean(mean(xcorrs(:,i,:,:,:),4),1)),[-.2,.2]) ; axis xy ; colormap jet; end
mxcorrs = squeeze(mean(mean(xcorrs,1),2)); 
allvox_mxcorrs = zeros(size(res_fmri,1),41); 
allvox_mxcorrs(maskinds,:) = mxcorrs; 
res_mxcorrs = reshape(allvox_mxcorrs,[size(mask.img,1),size(mask.img,2),size(mask.img,3),41]);
res_mxcorrs(isnan(res_mxcorrs)) = 0; 
cd E:\meanepis
ref41 = load_untouch_nii('ref41.nii.gz'); 
ref41.img = res_mxcorrs;
save_untouch_nii(ref41,'motion_xcorrs.nii.gz'); 
meanimg = load_untouch_nii('c:/shared/epireg/mean.nii.gz'); 
mrescorrs = double(squeeze(mean(res_mxcorrs(:,:,:,27:29),4))); 
for i=1:64
    subplottight(6,11,i);
    plotoverlayIntensity2D(squeeze(meanimg.img(i,:,:)),squeeze(mat2gray(abs(mrescorrs(i,:,:)))),(squeeze(mrescorrs(i,:,:))),90);
end
