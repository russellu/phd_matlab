clear all ; close all; 
subs = {'alex','dina','genevieve','jeremie'};
scans = {'den_retino_allstims_01','den_retino_allstims_02','den_retino_gamma_01','den_retino_gamma_02','den_retino_movie','den_retino_rest'};

sb = 1; 
cd(['c:/shared/coupling/',subs{sb}]); 

grid = load_untouch_nii('grid.nii.gz'); 
fmri = load_untouch_nii('grid_bp_clean_retino_movie.nii.gz');

cd(['c:/shared/coupling/',subs{sb},'/den_retino_movie']);
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

xcorrs = zeros(size(source_ts,1),50,41); 
for i=1:size(source_ts,1)
    for j=1:size(source_ts,2)
        source_ij = squeeze(source_ts(i,j,:)); 
        fmri_ij = squeeze(fmri_ts(i,:))'; 
        
        wsize = 50; wincr = 10; wcount=1;
        for w=100:wincr:length(source_ij)-100
            dyn_xcorrs(wcount,:) = xcorr(fmri_ij(w-wsize:w+wsize),source_ij(w-wsize:w+wsize),20,'coeff');
            wcount = wcount+1; 
        end
        xcorrs(i,j,:) = mean(dyn_xcorrs,1); 
    end
end
for i=1:200 ; subplottight(10,20,i);  imagesc(squeeze(xcorrs(i*10,:,:)),[-.2,.2]) ; set(gca,'XTick',[],'YTick',[]); axis xy;  end

ind = [8,31]; 
corrgrid = zeros(size(grid.img)); 
for i=1:length(gx)
    corrgrid(gx(i),gy(i),gz(i)) = xcorrs(i,ind(1),ind(2)); 
end
grid.img = corrgrid; 
cd(['c:/shared/coupling/',subs{sb}]); 
save_untouch_nii(grid,'corrgrid_movie.nii.gz'); 

%{
        % wsize = 50; wincr = 10; wcount=1;
        %for w=100:wincr:length(source_ij)-100
        %    scorrs(wcount) = corr2(source_ij(w-wsize:w+wsize),conved(w-wsize:w+wsize)); 
        %    wcount = wcount + 1; 
        %end
        %corrs(i,j) = mean(scorrs); 

%}


