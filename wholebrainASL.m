clear all ; close all ; 
cd c:/shared/asl2 ; ls 
subjs = {
    'nifti_gina_asl_dicom'   
    'nifti_vincent_asl_dicom'
    'nifti_russell2_asl_dicom' 
    'nifti_julie_asl_dicom' 
    'nifti_jeremie_asl_dicom'
    'nifti_alex_asl_dicom'    
} ;
stimtypes = [1,2,6] ; 
sub = 1 ; 
mask = load_untouch_nii('c:/shared/regasl/avg.nii.gz') ;
brainmask = mask.img > 2500 ; 
for sub=1:max(size(subjs)) ;
    cd(['c:/shared/asl2/',subjs{sub}]) ;
    asl=dir('al_reg_asl*') ; 
    bold=dir('al_reg_bold*') ; 
    
    % get the ASL
    disp(['loading asl images for ',subjs{sub}]) ;
    clear alld
    for i=1:max(size(asl))
        aslnii = load_untouch_nii(asl(i).name) ; 
        aslimg = aslnii.img ; 
        dasl = aslimg(:,:,:,1:60)-aslimg(:,:,:,61:120) ;
        alld(i,:,:,:,:) = dasl ; 
    end

    % get the BOLD
    disp(['loading bold images for ',subjs{sub}]) ;
    clear allb
    for i=1:max(size(bold))
        boldnii = load_untouch_nii(bold(i).name) ;  
        boldimg = boldnii.img ; 
        allb(i,:,:,:,:) = boldimg ; 
    end
    
    %%% get the triggers and canonical HRF convolved ideal time series
    disp(['calculating correlations for ',subjs{sub}]) ;
    boldtrigs = [2,4,6] ; asltrigs = [1,3,5] ; 
    bold_hrf = spm_hrf(2) ; % spm hrf with a tr of 2 seconds
    asl_hrf = spm_hrf(8) ; 
    stimTimes = dir('stimTimes*') ; clear times types stimvols ideal 
    for s=1:max(size(stimTimes)) ; 
        st = load(stimTimes(s).name) ;  
        st = st.stimTimes ; 
        for i=1:max(size(st)) ;
            times(s,i) = st{i}(1) ;
            types(s,i) = st{i}(2) ; 
        end
        boldtrvols = ceil(times(1,:)./2) ; allboldtrvols = ceil(times./2) ; 
        asltrvols = ceil(times(1,:)./8) ; allasltrvols = ceil(times./8) ; 
        
        % find the indices closest to the onset of the stimulus
        % ************ BOLD *************
        ideal = zeros(1,240) ; 
        for i=1:max(size(boldtrvols)) % fill the ideal with 1s where there was a stimulus
            ideal(boldtrvols(i):boldtrvols(i)+12) = 1 ; 
        end
        boldidealhrf = conv(ideal,bold_hrf) ; % convolve ideal with HRF    
        boldidealhrf = boldidealhrf(1:240) ; 
        
        
         % find the indices closest to the onset of the stimulus 
         % ********* ASL **********  
        ideal = zeros(1,60) ; 
        for i=1:max(size(asltrvols)) % fill the ideal with 1s where there was a stimulus
            ideal(asltrvols(i):asltrvols(i)+3) = 1 ; 
        end
        aslidealhrf = conv(ideal,asl_hrf) ; % convolve ideal with HRF    
        aslidealhrf = aslidealhrf(1:60) ;        
    end
    
    clear boldstimvols 
    stimcounts = [1,1,1] ;
    for i=1:max(size(boldtrigs)) 
        typesi = types(boldtrigs(i),:) ; 
        for j=1:max(size(typesi)) ;
            stimind = find(stimtypes==typesi(j)) ; 
            timeind = boldtrvols(j) ; 
            boldstimvols(i,stimind,stimcounts(stimind),:,:,:) = (squeeze(mean(allb(i,:,:,:,timeind+3:timeind+12),5))-squeeze(mean(allb(i,:,:,:,timeind-4:timeind),5))) ./ squeeze(mean(allb(i,:,:,:,timeind-4:timeind),5)) ; 
            stimcounts(stimind) = stimcounts(stimind) + 1 ; 
        end
    end
    mbstv = squeeze(mean(mean(boldstimvols,1),3)) ; 
    mbstv(isnan(mbstv)) = 0 ; mbstv(isinf(mbstv)) = 0 ;
    for i=1:3
       mbstv(i,:,:,:) = squeeze(mbstv(i,:,:,:)).*brainmask ;  
       mbstv(i,:,:,:) = medfilt3(squeeze(mbstv(i,:,:,:))) ; 
    end
    imagesc(squeeze(mbstv(1,:,:,18)),[-.025,.025])   

    clear aslstimvols 
    stimcounts = [1,1,1] ;
    for i=1:max(size(asltrigs)) 
        typesi = types(asltrigs(i),:) ; 
        for j=1:max(size(typesi)) ;
            stimind = find(stimtypes==typesi(j)) ; 
            timeind = asltrvols(j) ; 
            aslstimvols(i,stimind,stimcounts(stimind),:,:,:) = squeeze(mean(alld(i,:,:,:,timeind+1),5)) ./ squeeze(mean(alld(i,:,:,:,timeind-1),5)) ; 
            stimcounts(stimind) = stimcounts(stimind) + 1 ; 
        end
    end
    mdstv = squeeze(mean(mean(aslstimvols,1),3)) ; 
    mdstv(isnan(mdstv)) = 0 ; mdstv(isinf(mdstv)) = 0 ;
    for i=1:3
       mdstv(i,:,:,:) = squeeze(mdstv(i,:,:,:)).*brainmask ;  
       mdstv(i,:,:,:) = medfilt3(squeeze(mdstv(i,:,:,:))) ; 
    end   
    %imagesc(squeeze(mdstv(1,:,:,15)),[-1,1])   
    alpha = 0.6 ;  
    beta = 1.3 ; 
	M=0.11; 
    cmro2t = ((1-mbstv./M).^(1/beta)).*(mdstv.^(1-alpha/beta)) ; 
    cmro2t = real(cmro2t) ; cmro2t(isnan(cmro2t)) = 0 ; cmro2t(isinf(cmro2t)) = 0 ; 
    %for i=1:10 ; subplot(3,4,i) ; imagesc(squeeze(cmro2t(3,:,:,10+i))) ;  end
    allcmro2t(sub,:,:,:,:) = cmro2t ; 
end
for i=1:10 ; subplot(3,4,i) ; imagesc(squeeze(mean(allcmro2t(:,3,:,:,10+i),1))) ; end


