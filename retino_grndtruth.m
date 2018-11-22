% estimate hrf from event related designclear all ; close all ; 
% stimulus starts in the bottom left quadrant (225deg=S1). 
clear all ; close all ;
trigs{1} = {'S 11','S 12','S 13','S 14'} ; 
trigs{2} = {'S 21','S 22','S 23','S 24'} ; 
trigs{3} = {'S 31','S 32','S 33','S 34'} ; 
trigs{4} = {'S 41','S 42','S 43','S 44'} ; 
trigs{5} = {'S 51','S 52','S 53','S 54'} ; 
trigs{6} = {'S 61','S 62','S 63','S 64'} ; 
subs= {'alex','dina','genevieve','jeremie','karl','russell','valerie'} ;

for subj=6%:length(subs)
    sub = subs{subj} ; 
    cd(['c:/shared/badger_mri/',sub,'/nii']) ; ls 
    fmris=dir('bp_reg_mc*allstim*nii*') ; 
    cd(['c:/shared/badger_eeg/',sub]) ; 
    eegs = dir('1hz*allstim*set') ; 

    stimcounts = [1,1,1] ;  
    for f=1:length(fmris) ; 
       cd(['c:/shared/badger_mri/',sub,'/nii']) ; ls 
       fmri = load_untouch_nii(fmris(f).name) ;  
       fimg = fmri.img ; 
       if f==1 ; allstims = zeros(6,size(fimg,1),size(fimg,2),size(fimg,3),round(20./0.693)) ; end

       TR = fmri.hdr.dime.pixdim(5) ; 
       cd(['c:/shared/badger_eeg/',sub]) ; ls 
       eeg = pop_loadset(eegs(f).name) ; 
       events = {eeg.urevent.type} ; 
       onsetinds = find(strcmp('R128',events)) ; 
       lats = {eeg.urevent.latency} ; 
       onsetlats = cell2mat(lats(onsetinds)) ; 
       vlats = cell2mat(lats(onsetinds)) ; 
       trtopeak = find(spm_hrf(TR)==max(spm_hrf(TR))) ; 
       for i=1:length(trigs) ; % for all orientations
              for j=1:length(trigs{i})      
                  trigi = find(strcmp(trigs{i}(j),events)) ; 
                  latij = cell2mat(lats(trigi)) ; 
                  diffs = onsetlats-latij ; 
                  indij = find(abs(diffs) == min(abs(diffs))) ; 
                  allstims(i,:,:,:,:) = squeeze(allstims(i,:,:,:,:)) + fimg(:,:,:,indij:indij+round(20/TR)-1) ; 
              end
       end
    end

    hrf = spm_hrf(0.693) ; 
    hrf = hrf(1:round(10/TR)) ; 

    % shift the HRF and correlate
    stimperiod = zeros(1,round(12/TR)) ; 
    clear spikes convspikes
    for i=1:length(stimperiod)
        spikes(i,:) = stimperiod ; 
        spikes(i,i) = 1 ;       
        convspikes(i,:) = conv(spikes(i,:),hrf) ; 
    end
    convspikes = squeeze(convspikes(:,1:size(allstims,5))) ; 
    clear voxcorrs ; 
    for i=1:6
        for j=1:size(convspikes,1)
            voxcorrs(i,:,:,:,j) = voxcorr(squeeze(allstims(i,:,:,:,:)),squeeze(convspikes(j,:))) ; 
        end
    end
    
    as(:,:,:,1) = squeeze(voxcorrs(5,:,:,:,1)) ; 
    as(:,:,:,2) = squeeze(voxcorrs(4,:,:,:,1)) ; 
    as(:,:,:,3) = squeeze(voxcorrs(3,:,:,:,1)) ; 
    as(:,:,:,4) = squeeze(voxcorrs(2,:,:,:,1)) ; 
    as(:,:,:,5) = squeeze(voxcorrs(1,:,:,:,1)) ; 
    as(:,:,:,6) = squeeze(voxcorrs(6,:,:,:,1)) ; 
    
    ref = load_untouch_nii('allstim_ref6.nii.gz') ; 
    ref.img = as ; mask = load_untouch_nii('anatpos.nii.gz')  ;
    for i=1:size(ref.img,4) ; ref.img(:,:,:,i) = ref.img(:,:,:,i).*mask.img ; end
    save_untouch_nii(ref,'as6.nii.gz') ; 
    cd(['c:/shared/badger_mri/',sub,'/nii']) ; ls 
    allstimref = load_untouch_nii('allstim_ref.nii.gz') ; 
    for i=1:6
       % ref = load_untouch_nii('f_mc_retino_allstims_01.nii.gz') ; 
       % corrs = voxcorr(squeeze(allstims(i,:,:,:,:)),hrf') ;              
       % ref.img = corrs ; 
       % save_untouch_nii(ref,[num2str(i),'_allstims.nii.gz']) ; 
       allstimref.img = squeeze(voxcorrs(i,:,:,:,:)) ; 
       save_untouch_nii(allstimref,['tcourse_',num2str(i),'.nii.gz']) ; 
    end
    
    
end
