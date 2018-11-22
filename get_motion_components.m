clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 
source_folders = {'den_retino_allstims_01','den_retino_allstims_02','den_retino_gamma_01','den_retino_gamma_02','den_retino_movie','den_retino_rest'};  
clean_fmri_names = {'bp_clean_retino_allstims_01','bp_clean_retino_allstims_02','bp_clean_retino_gamma_01','bp_clean_retino_gamma_02','bp_clean_retino_movie','bp_clean_retino_rest'};  
gsr_fmri_names = {'bp_gsr_clean_retino_allstims_01','bp_gsr_clean_retino_allstims_02','bp_gsr_clean_retino_gamma_01','bp_gsr_clean_retino_gamma_02','bp_gsr_clean_retino_movie','bp_gsr_clean_retino_rest'}; 
eeg_set_names = {'retino_allstims_01','retino_allstims_02','retino_gamma_01','retino_gamma_02','retino_movie','retino_rest'}; 
scans = {'mc_retino_allstims_01','mc_retino_allstims_02','mc_retino_gamma_01','mc_retino_gamma_02','mc_retino_movie','mc_retino_rest'}; 
mask = load_untouch_nii('c:/shared/epireg/meanmask.nii.gz'); 
mean_nii = load_untouch_nii('c:/shared/epireg/mean.nii.gz'); 

maskinds = find(mask.img==1); 
for sb=1:length(subs)    
    corrs = load_untouch_nii(['e:/fmris/badger_',subs{sb},'/atlas_gamma_mcorrs.nii.gz']);
    corrs.img(isnan(corrs.img)) = 0; 
    figure,
    for srcf=1:6   
        cd(['E:\rawbadger\badger_mri\',subs{sb},'\nii\',scans{srcf},'.nii.mat']);
        mats = dir('MAT*');
        ls 
        
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

        cd(['E:/badger_eeg/',subs{sb}]);
        eeg = pop_loadset([eeg_set_names{srcf},'.set']);    
        
        trigs = {eeg.urevent.type};
        trtrigs = find(strcmpi('R128',trigs)); 
        lats = {eeg.urevent.latency}; 
        trlats = cell2mat(lats(trtrigs)); 
        
        allweights = load('allweights'); allweights = allweights.allweights;        
        acts = allweights*eeg.data; 
        
        freqs = 1:10:130; 
        clear res_highf
        for f=1:length(freqs)
            highf = eegfiltfft(acts,eeg.srate,freqs(f)-2,freqs(f)+2); 
            abs_highf = abs(highf); 
            abs_highf = abs_highf(:,trlats(1):trlats(end)); 
            res_highf(f,:,:) = imresize(abs_highf,[64,length(trlats)]); 
        end
        
        allweights_sortmotion = load('allweights_sortmotion');
        allweights_sortmotion = allweights_sortmotion.allweights_sortmotion; 
        subplot(2,3,srcf);
        
        comps = allweights_sortmotion(1:5); 
        clear z_filt_highf
        for i=1:5
           z_filt_highf(:,i,:) = zscore(eegfiltfft(squeeze(res_highf(:,comps(i),:)),1/0.693,0.005,1),[],2); 
        end

        zs = zscore(squeeze(mean(mean(z_filt_highf(8:end,:,:),2),1)));
        [sv,si]  =sort(zs,'descend'); 
        plot(zs) ; vline(si(1:round(length(si)*0.1))); 
        
        save(['si_',eeg_set_names{srcf}],'si')
        
        for i=1:64
            for f=1:size(res_highf,1)
                xcorrs(sb,srcf,i,f,:) = xcorr(zsummat(20:size(res_highf,3)-21),squeeze(res_highf(f,i,21:end-20)),20,'coeff') ;
            end  
        end      
  
    end
    
    %[sv,si] = sort(squeeze(mean(mean(mean(xcorrs(sb,:,:,10:end,20:22),4),5),2)),'descend'); 
    %subplot(2,4,sb); imagesc(squeeze(mean(mean(xcorrs(sb,:,si(1:3),:,:),2),3)),[-.2,.2]); colormap jet;    
    %allweights_sortmotion = si; save('allweights_sortmotion','allweights_sortmotion'); 
end



