clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','tegan','valerie'}; 
scans = {'*allstim*01*','*allstim*02*','*gamma*01*','*gamma*02*','*movie*','*rest*'};
freqs = {'delta','theta','alpha','beta','mid_gamma','high_gamma'};
rsns = {'anterior_visual.nii','dmn.nii','higher_visual.nii','primary_visual.nii','motor.nii'};

cd C:\shared\epireg\rsns\networks;

for r=1:length(rsns)
   rsn = load_untouch_nii(rsns{r}); 
   [sv,si] = sort(rsn.img(:),'descend'); 
   rsnvox(r,:) = si(1:1000);    
end

for sb=1:length(subs) 
    cd(['e:/sim/fmri/',subs{sb}]);
    meanimg = load_untouch_nii('mean.nii.gz'); 
    mask = meanimg.img > mean(meanimg.img(:)*1.2); 
    maskvox = find(mask == 1); 
    allmots = load('allmots'); allmots = allmots.allmots; 
    
    for sc=1:length(scans)
        fmri = dir(['bp_gsr_clean_',scans{sc},'.nii']);
        fmri = load_untouch_nii(fmri.name); 
        resfmri = reshape(fmri.img,[numel(fmri.img(:,:,:,1)),size(fmri.img,4)]); 
        resfmrivox = resfmri(maskvox,:); 
        
        motion_i = allmots{sc}; 
        [sv,si] = sort(motion_i,'descend'); 
        
        badinds = si(1:length(si)/10); 
        goodinds = si(length(si)/10:end); 
        hrf = spm_hrf(0.693); 
        for fr=1:length(freqs)
            eeg = dir(['atlas_',scans{sc},freqs{fr},'*nii']);
            eeg = load_untouch_nii(eeg.name); 
            reseeg = reshape(eeg.img,[numel(eeg.img(:,:,:,1)),size(eeg.img,4)]); 
            reseeg(:,badinds) = repmat(mean(reseeg(:,goodinds),2),[1,length(badinds)]); 
            clear eegrsns fmrirsns
            for i=1:length(rsns)
                eegrsn = smooth(eegfiltfft(squeeze(mean(reseeg(rsnvox(i,:),:),1)),1/0.693,0.001,2))'; 
                %eegrsn = smooth((squeeze(mean(reseeg(rsnvox(i,:),:),1))),30)';
                %eegrsn = eegrsn(1:size(reseeg,2)); 
                fmrirsn = squeeze((resfmri(rsnvox(i,:),:)));                 
                %xcorrs(sb,sc,fr,i,:) = xcorr(fmrirsns(i,35:size(eegrsns,2)-35),eegrsns(i,35:end-35),25,'coeff');           
                dyncorrs = zeros(1,size(fmrirsn,2)); 
                for j=45:size(fmrirsn,2)-45
                    dyncorrs(j) = mean(mean(corr(fmrirsn(:,j-15:j+15)',fmrirsn(:,j-15:j+15)')));   
                end
                %xcorrs(sb,sc,fr,i,:) = xcorr(dyncorrs(45:size(eegrsn,2)-45),eegrsn(45:end-45),25,'coeff'); 
                
                dyncorrs = dyncorrs(45:size(eegrsn,2)-45); 
                eegrsn = eegrsn(45:end-45); 
                
                [sv,si] = sort(dyncorrs,'descend'); 
                sortpower(sb,sc,fr,i,1:200) = (eegrsn(si(1:200))); 
                sortpower(sb,sc,fr,i,201:400) = (eegrsn(si(end-199:end)));     
                
                straightcorrs(sb,sc,fr,i) = corr(eegrsn',dyncorrs'); 
            end
            disp([subs{sb},' ',freqs{fr},' ',scans{sc}]);
        end   
    end

end
straightcorrs(isnan(straightcorrs)) = 0; 
cd E:\sim\saved
%mc_yes_gsr = freqzimgs ; save('mc_yes_gsr','mc_yes_gsr'); 
%mc_no_gsr_xcorrs = xcorrs; save('mc_no_gsr_xcorrs','mc_no_gsr_xcorrs');
%rsn_corr_nogsr = xcorrs; save('rsn_corr_nogsr','rsn_corr_nogsr'); 
sortpower_gsr = sortpower; save('sortpower_gsr','sortpower_gsr'); 
straightcorrs_gsr = straightcorrs ; save('straightcorrs_gsr','straightcorrs_gsr'); 




