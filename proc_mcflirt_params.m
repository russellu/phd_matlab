clear all ; close all ; 

subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 
scans = {'mc_retino_allstims_01','mc_retino_allstims_02','mc_retino_gamma_01','mc_retino_gamma_02','mc_retino_movie','mc_retino_rest'}; 
eegscans = {'retino*allstim*01*set','retino*allstim*02*set','retino*gamma*01*set','retino*gamma*02*set','retino*movie*set','retino*rest*set'}; 

for sb=8%:length(subs)
    cd(['E:\rawbadger\badger_mri\',subs{sb},'\nii']);   
    for sc=6%:length(scans)
        cd(['E:\rawbadger\badger_mri\',subs{sb},'\nii\',scans{sc},'.nii.mat']);
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
        
        cd(['E:\badger_eeg\',subs{sb}]);  
        eegscan = dir(eegscans{sc}); 
        disp(eegscan.name); 
        
        eeg = pop_loadset(eegscan.name); 
        
        eeg.data = eegfiltfft(eeg.data,eeg.srate,1,128); 
        [weights,sphere] = runica(eeg.data,'maxsteps',128); 
        eeg.data = weights*sphere*eeg.data; 
           
        raweeg = eeg; 
        trigs = {eeg.urevent.type};
        trtrigs = find(strcmpi('R128',trigs)); 
        lats = {eeg.urevent.latency}; 
        trlats = cell2mat(lats(trtrigs)); 
          
        freqs = 1:2:100; 
        clear res_filts bcg_res_filts    
        for f=1:length(freqs)
           filt_f = eegfiltfft(eeg.data,eeg.srate,freqs(f)-1,freqs(f)+1);  
           filt_f = filt_f(:,trlats(1):trlats(end)); 
           res_filts(f,:,:) = imresize(abs(filt_f),[64,length(trlats)]);  
        end
        
        res_filts(:,:,1) = []; 
        for i=1:50
           res_filts(i,:,:) = eegfiltfft(squeeze(res_filts(i,:,:)),1/0.693,0.005,2);  
        end
        zsummat = zsummat(1:size(res_filts,3)); 
        
        for i=1:50
            for j=1:64
                xcorrs(sb,sc,i,j,:) = xcorr(zsummat(20:end-20),squeeze(res_filts(i,j,20:end-20)),20,'coeff'); 
            end
        end        
    end
end


%for i=1:6 ; subplot(2,3,i) ; imagesc(squeeze(mean(mean(xcorrs(:,i,:,:,:),4),1)),[-.2,.2]) ; axis xy ; colormap jet; end

