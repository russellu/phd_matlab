clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','tegan','valerie'}; 
scans = {'*allstim*01*','*allstim*02*','*gamma*01*','*gamma*02*','*movie*','*rest*'};
freqs = {'delta','theta','alpha','beta','mid_gamma','high_gamma'};

for sb=1:length(subs) 
    cd(['e:/sim/fmri/',subs{sb}]);
    meanimg = load_untouch_nii('mean.nii.gz'); 
    mask = meanimg.img > mean(meanimg.img(:)*1.2); 
    maskvox = find(mask == 1); 
    allmots = load('allmots'); allmots = allmots.allmots; 
    
    for sc=1:length(scans)
        fmri = dir(['bp_clean_',scans{sc},'.nii']);
        fmri = load_untouch_nii(fmri.name); 
        resfmri = reshape(fmri.img,[numel(fmri.img(:,:,:,1)),size(fmri.img,4)]); 
        resfmrivox = resfmri(maskvox,:); 
        
        motion_i = allmots{sc}; 
        %subplot(2,1,1); plot(motion_i); hold on; 
        [sv,si] = sort(motion_i,'descend'); 
        
        badinds = si(1:length(si)/10); 
        goodinds = si(length(si)/10:end); 
        
        %motion_i(badinds) = mean(motion_i(goodinds)); 
         %plot(motion_i); xlabel('time (TR=0.693s)'); ylabel('high gamma (100+Hz) amp'); 
         %legend({'before interp.','after interp.'});
        
        
        for fr=1:length(freqs)
            eeg = dir(['atlas_',scans{sc},freqs{fr},'*nii']);
            eeg = load_untouch_nii(eeg.name); 
            reseeg = reshape(eeg.img,[numel(eeg.img(:,:,:,1)),size(eeg.img,4)]); 
            reseegvox = reseeg(maskvox,:); 
            reseegvox(:,badinds) = repmat(mean(reseegvox(:,goodinds),2),[1,length(badinds)]); 
            
            hrf = spm_hrf(0.693); 
            corrs = zeros(1,size(resfmrivox,1)); 
            for v=1:size(reseegvox,1)
                conved_v = conv(reseegvox(v,:),hrf,'full');
                conved_v = conved_v(1:size(reseegvox,2)); 
                corrs(v) = corr(conved_v(30:size(reseegvox,2)-30)',resfmrivox(v,30:size(reseegvox,2)-30)'); 
            end
            zimg = zeros(size(meanimg.img));
            zimg(maskvox) = corrs; 
            zimg(isnan(zimg)) = 0; 
            
            freqzimgs(:,:,:,fr,sc,sb) = zimg; 
            disp([subs{sb},' ',scans{sc}]);
            
        end   
    end
end

cd E:\sim\saved
mc_no_gsr = freqzimgs ; save('mc_no_gsr','mc_no_gsr'); 






