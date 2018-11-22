clear all ; close all ; 
stimtypes = [1,2,6] ; 

basetime = 2 ; 
tasktime = 14 ; 

convbasetime = 0 ; 
convtasktime = 16 ; 

sub = 1 ; 
    cd('c:/shared/short_tr/mb_russ2_img') ;
    %asl=dir('reg_asl*') ; 
   %
    bold=dir('reg*') ; boldtrigs = [1,3] ; boldTR = 0.41 ; stitle = 'multi-band (TR=0.41s), #trials = 20' ;
    %bold=dir('mc*TR275*') ; boldtrigs = [2,4] ; boldTR = 0.275 ; stitle = 'single-band few slices (TR=0.41s), #trials = 20' ;   
    % get the BOLD
    clear allb
    for i=1:max(size(bold))
        boldnii = load_untouch_nii(bold(i).name) ;  
        boldimg = boldnii.img ; 
        allb(i,:,:,:,:) = boldimg ; 
    end
    bold_hrf = spm_hrf(boldTR) ; % spm hrf with a tr of 2 seconds
    bold_hrf(2:round(2/boldTR)) = -.03 ;
    bold_hrf = smooth(bold_hrf,10) ; 
    stimTimes = dir('stimTimes*') ; clear times types stimvols ideal 
    for s=1:max(size(stimTimes)) ; 
        st = load(stimTimes(s).name) ;  
        st = st.stimTimes ; 
        for i=1:max(size(st)) ;
            times(s,i) = st{i}(1) ;
            types(s,i) = st{i}(2) ; 
        end
        boldtrvols = round(times(1,:)./boldTR) ; allboldtrvols = round(times./boldTR) ; 
        
        % find the indices closest to the onset of the stimulus
        % ************ BOLD *************
        ideal = zeros(1,size(allb,5)) ; 
        for i=1:max(size(boldtrvols)) % fill the ideal with 1s where there was a stimulus
            ideal(boldtrvols(i):boldtrvols(i)+round(2/boldTR)) = 1 ; 
        end
        boldidealhrf = conv(ideal,bold_hrf) ; % convolve ideal with HRF    
        boldidealhrf = boldidealhrf(1:size(allb,5)) ; 
    end 
    %%%% epoch the raw BOLD time series 
    clear boldepochs fullboldepochs
    for sc=1:size(allb,1) %  for all scans
        boldtrigvols = allboldtrvols(boldtrigs(sc),:) ; % stimulus onset for that scan
        boldtrigtypes = types(boldtrigs(sc),:) ; 
        for i=1:max(size(stimtypes))
           boldtypinds = boldtrigvols(find(boldtrigtypes==stimtypes(i))) ;  
           for j=1:max(size(boldtypinds))
                fullboldepochs(i,(sc*10-10)+j,:,:,:,:) = allb(sc,:,:,:,boldtypinds(j)-round(basetime./boldTR):boldtypinds(j)+round(tasktime./boldTR)) ;         
                boldepochs(i,(sc*10-10)+j,:,:,:,:) = allb(sc,:,:,:,boldtypinds(j)-round(convbasetime./boldTR):boldtypinds(j)+round(convtasktime./boldTR)) ;                
           end
        end
    end
    
    refvol = load_untouch_nii('mb1.nii.gz') ; 

    % *********** bold corrvols *************
    boldcorrvol = zeros(1,1,1,size(boldepochs,6)) ; boldcorrvol(1,1,1,:) = bold_hrf(1:size(boldepochs,6)) ; 
    boldcorrvol = repmat(boldcorrvol,[size(allb,2),size(allb,3),size(allb,4),1]) ; 
    for i=1:size(boldepochs,1) ; disp(i) ; 
        for j=1:size(boldepochs,2)
            allcorrs(i,j,:,:,:) = voxcorr(boldcorrvol,squeeze(boldepochs(i,j,:,:,:,:))) ;          
        end
    end
    refvol.img = squeeze(mean(mean(allcorrs))) ; 
    save_untouch_nii(refvol,'dipcorrs.nii.gz') ; 
    
    mtrials = squeeze(mean(allcorrs,2)) ; 
    mmtrials = squeeze(mean(mtrials)) ; 
    
    corrthresh = 0.2 
    inds = find(mmtrials>corrthresh) ; 
    voxcount = 1 ; clear corrvox ;
    % get the correlating voxels
    for i=1:size(boldepochs,3) 
        for j=1:size(boldepochs,4)
            for k=1:size(boldepochs,5)
                if mmtrials(i,j,k) > corrthresh
                    corrvox(voxcount,:,:,:) = squeeze(fullboldepochs(:,:,i,j,k,:)) ; 
                    voxcount = voxcount + 1 ; 
                end
            end
        end
    end
    % baseline correct the correlating voxels
    bases = zeros(size(corrvox,1),size(corrvox,2),size(corrvox,3),1) ; 
    bases(:,:,:,1) = squeeze(mean(corrvox(:,:,:,1:round(basetime./boldTR)),4)) ; bases = repmat(bases,[1,1,1,size(corrvox,4)]) ;
    bcorrvox = corrvox - bases ; mbcorrvox = squeeze(mean(bcorrvox,1)) ; 
    figure,
    dipt = round(basetime/boldTR:basetime/boldTR+2/boldTR) ; 
    subplot(2,2,1) ; errorbar(squeeze(mean(mbcorrvox([1,2],:,:),2))',squeeze(std(mbcorrvox([1,2],:,:),0,2))'./sqrt(20),'LineWidth',2) ;  hline(0,'k') ; %vline(round(1./boldTR),'r') ; vline(round(1./boldTR)+round(2./boldTR),'r') ; 
    vline(dipt(1),'r'),vline(dipt(max(size(dipt))),'r') ;     
    subplot(2,2,2) ; errorbar(squeeze(mean(mbcorrvox([1,2],:,1:max(dipt)+round(2/boldTR)),2))',squeeze(std(mbcorrvox([1,2],:,1:max(dipt)+round(2/boldTR)),0,2))'./sqrt(20),'LineWidth',3) ;  hline(0,'k') ; %vline(round(1./boldTR),'r') ; vline(round(1./boldTR)+round(2./boldTR),'r') ; 
    vline(dipt(1),'r'),vline(dipt(max(size(dipt))),'r') ;  
    
    dip = squeeze(mean(mbcorrvox(:,:,dipt),3)) ; %legend({'unperturbed','5%contrast','60%random'}) ;
    
    
    
    
    
    
 
   % subplot(2,2,2) ; barwitherr(std(dip,0,2)./sqrt(20),mean(dip,2)) ; set(gca,'XTickLabel',{'unperturbed','5%contrast','60%random'}) ; suptitle('strength of initial dip') ;  
    
    
    %{
    rgb(:,:,:,1) = squeeze(mean(mean(allcorrs))) ; 
    rgb(:,:,:,3) = fullim ; 
    for i=4:16 ; subplot(3,4,i-3) ; 
        imagesc(squeeze(uint8(mat2gray(squeeze(rgb(:,:,i,:)))*255))) ; 
    end
    %}
    
    
    
    
    
    
    
    
    
    
    
    