clear all ; close all; 
subs = {'alex','dina','genevieve','jeremie'};
scans = {'den_retino_allstims_01','den_retino_allstims_02','den_retino_gamma_01','den_retino_gamma_02','den_retino_movie','den_retino_rest'};
fscans = {'grid_bp_clean_retino_allstims_01.nii.gz','grid_bp_clean_retino_allstims_02.nii.gz','grid_bp_clean_retino_gamma_01.nii.gz'...
    'grid_bp_clean_retino_gamma_02.nii.gz','grid_bp_clean_retino_movie.nii.gz','grid_bp_clean_retino_rest.nii.gz'}; 

sb = 3; 
for scn=1:6
cd(['c:/shared/coupling/',subs{sb}]); 

grid = load_untouch_nii('grid.nii.gz'); 
fmri = load_untouch_nii([fscans{scn}]);

cd(['c:/shared/coupling/',subs{sb},'/',scans{scn}]);
hzs = dir('*gz'); 
hz1 = load_untouch_nii(hzs(1).name); 
hzimgs = zeros(size(grid.img,1),size(grid.img,2),size(grid.img,3),50,size(hz1.img,4)); 
for hz=1:length(hzs)
    hzi = load_untouch_nii(hzs(hz).name); 
    hzimgs(:,:,:,hz,:) = hzi.img; 
end

gridinds = find(grid.img==1) ;
[gx,gy,gz] = ind2sub(size(grid.img),gridinds); 
% create vectors of FMRI signal and EEG signal

fmri_ts = zeros(length(gx),size(hzimgs,5)); 
source_ts = zeros(length(gx),size(hzimgs,4),size(hzimgs,5)); 

for i=1:length(gx)
    fmri_ts(i,:) = squeeze(fmri.img(gx(i),gy(i),gz(i),1:size(hzimgs,5))); 
end
for i=1:length(gx)
    source_ts(i,:,:) = squeeze(hzimgs(gx(i),gy(i),gz(i),:,:)); 
end
hrf = spm_hrf(0.693); 
corrs = zeros(size(source_ts,1),size(source_ts,2)); 
for i=1:size(source_ts,1)
    for j=1:size(source_ts,2)
        fmri_i = fmri_ts(i,1:size(source_ts,3))';
        source_ij = squeeze(source_ts(i,j,:)); 
        conved = conv(source_ij,hrf,'full');
        conved = conved(1:size(source_ts,3)); 
        
        wsize = 50; wincr = 10; wcount=1;
        for w=100:wincr:length(source_ij)-100
            scorrs(wcount) = corr2(fmri_i(w-wsize:w+wsize),conved(w-wsize:w+wsize)); 
            wcount = wcount + 1; 
        end
        corrs(i,j) = mean(scorrs); 
        
        
        %corrs(i,j) = corr2(conved,fmri_i); 
    end
end

corrgrid = zeros(size(grid.img,1),size(grid.img,2),size(grid.img,3),50); 
for i=1:length(gx)
    corrgrid(gx(i),gy(i),gz(i),:) = corrs(i,:); 
end
cd(['c:/shared/coupling/',subs{sb}]); 
gridnii = make_nii(corrgrid); gridnii.hdr.pixdim(2:4) = 8; 
save_nii(gridnii,['hrfcorrs_',scans{scn},'.nii.gz']); 
end


