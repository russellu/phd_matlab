clear all ; close all ; 

subs = {'alex','dina','genevieve','jeremie','russell','tegan','valerie'}; 
scans = {'mc_retino_allstims_01','mc_retino_allstims_02','mc_retino_gamma_01','mc_retino_gamma_02','mc_retino_movie','mc_retino_rest'}; 
eegscans = {'retino*allstim*01*set','retino*allstim*02*set','retino*gamma*01*set','retino*gamma*02*set','retino*movie*set','retino*rest*set'}; 
fmris = {'bp_clean_retino_allstims_01.nii.gz','bp_clean_retino_allstims_02.nii.gz','bp_clean_retino_gamma_01.nii.gz',...
    'bp_clean_retino_gamma_02.nii.gz','bp_clean_retino_movie.nii.gz','bp_clean_retino_rest.nii.gz'};
mask = load_untouch_nii('c:/shared/epireg/meanmask.nii.gz'); 
maskinds = find(mask.img==1); 

xcorrs = zeros(8,6,35994,41); 
for sb=1:length(subs)
    for sc=1:length(scans)
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
        
        
        allzsums{sc} = zsummat; 
        
        %{
        zsummat = smooth(zsummat) ;
        zsummat = eegfiltfft(zsummat',1/0.693,0.005,1); 
        hrf = spm_hrf(0.693); 
        conv_motion = conv(zsummat,hrf,'full'); 
        conv_motion = conv_motion(1:length(zsummat)); 
        %}

        % get the BOLD
        %{
       cd(['E:\fmris\badger_',subs{sb},'\atlas_fmri']);
       bold = load_untouch_nii(fmris{sc}); 
       
       hrfcorrs(:,:,:,sb,sc) = voxcorr(bold.img(:,:,:,30:end-31),conv_motion(30:end-30)); 
       nohrfcorrs(:,:,:,sb,sc) = voxcorr(bold.img(:,:,:,30:end-31),zsummat(30:end-30)); 
        %}
    end
    cd(['E:\sim\fmri\',subs{sb}]);
    save('allzsums','allzsums'); 
end



