clear all ; close all; 
subs = {'alex','dina','genevieve','jeremie','russell','tegan','valerie'}; 
sets = {'retino_allstim*01*set','retino_allstim*02*set','retino_gamma*01*set','retino_gamma*02*set','retino_movie*set','retino_rest*set'}; 
fmris = {'bp_clean_retino_allstims_01.nii.gz','bp_clean_retino_allstims_02.nii.gz','bp_clean_retino_gamma_01.nii.gz',...
    'bp_clean_retino_gamma_02.nii.gz','bp_clean_retino_movie.nii.gz','bp_clean_retino_rest.nii.gz'};

fmris = {'bp_gsr_clean_retino_allstims_01.nii.gz','bp_gsr_clean_retino_allstims_02.nii.gz','bp_gsr_clean_retino_gamma_01.nii.gz',...
    'bp_gsr_clean_retino_gamma_02.nii.gz','bp_gsr_clean_retino_movie.nii.gz','bp_gsr_clean_retino_rest.nii.gz'};

freqs = {[1,3],[5,15],[16,25],[35,60]}; 

for sb=1:length(subs)
        
    % get the bcg (and electrode bcg) 
    for st=6%:length(sets)
       cd(['E:\badger_eeg\',subs{sb}]);
       set_st = dir(sets{st}); 
       eeg = pop_loadset(set_st.name); 
       
       trigs = {eeg.urevent.type};
       r128s = find(strcmpi(trigs,'R128')); 
       lats = cell2mat({eeg.urevent.latency}); 
       r128lats = lats(r128s);     
       % filter the bcg and downsample to number of BOLD samples 
       
       bcgica = load('bcgica'); bcgica = bcgica.bcgica; 
       weights = bcgica{1}; sphere = bcgica{2}; 
       acts = weights*sphere*eeg.data; 
       %eeg.data = acts; 
       
       filts = zeros(length(freqs),64,size(eeg.data,2)); 
       for f=1:length(freqs)
           filts(f,:,:) = abs(eegfiltfft(eeg.data,eeg.srate,freqs{f}(1),freqs{f}(2))); 
       end
       filts = filts(:,:,r128lats(1):r128lats(end)); 
       clear resfilts convedfilts
       for i=1:length(freqs)
          resfilts(i,:,:) = imresize(squeeze(filts(i,:,:)),[64,length(r128lats)]);  
       end
       hrf = spm_hrf(0.693); 
       for i=1:length(freqs)
          for j=1:64
             conved = conv(squeeze(resfilts(i,j,:)),hrf,'full'); 
             conved = conved(1:size(resfilts,3)); 
             convedfilts(i,j,:) = conved; 
          end
       end
        
       % get the BOLD
       cd(['E:\fmris\badger_',subs{sb},'\atlas_fmri']);
       bold = load_untouch_nii(fmris{st}); 
       
       for i=1:length(freqs)
          boldcorrs(:,:,:,i,sb) = voxcorr(bold.img(:,:,:,30:length(r128lats)-30),squeeze(convedfilts(i,2,30:length(r128lats)-30))); 
       end
       
    end

end

figure,
for i=1:4 ; subplot(2,3,i)
   imagesc(squeeze(mean(mean(boldcorrs(:,:,12:16,i,:),5),3)),[-.1,.1]);  
end

figure,
cd E:\sim\fmri\alex ; meanatlas = load_untouch_nii('mean.nii.gz'); 
slices = [5,8,10,12,14,16,20]; 
for sl = 1:length(slices)
   subplottight(1,7,sl); 
   corrimg = squeeze(mean(boldcorrs(:,:,slices(sl),2,:),5)); corrimg(1,1) = -.2; corrimg(1,2) = .2; corrimg(corrimg>.2) = 0.2; corrimg(corrimg<-.2) = -.2; 
   corrimg(isnan(corrimg)) = 0; 
   plotoverlayIntensity2D(meanatlas.img(:,:,slices(sl)),mat2gray(abs(corrimg)),(corrimg),270);
    
end








