clear all ; close all ; 
stimtypes = [1,2,6] ; 
sub = 1 ; 
    cd(['c:/shared/short_tr/mb_russ2_img']) ;
    %asl=dir('reg_asl*') ; 
   %
  bold=dir('mc*TR410*') ; boldtrigs = [1,3] ; boldTR = 0.41 ; stitle = 'multi-band (TR=0.41s), #trials = 20' ;
   %bold=dir('*TR275*') ; boldtrigs = [2,4] ; boldTR = 0.275 ; stitle = 'single-band few slices (TR=0.41s), #trials = 20' ;

    
    % get the BOLD
    clear allb
    for i=1:max(size(bold))
        boldnii = load_untouch_nii(bold(i).name) ;  
        boldimg = boldnii.img ; 
        allb(i,:,:,:,:) = boldimg ; 
    end
    %size(allb) 
    %disp('filtering data') ; 
    %resb = reshape(allb,[2*64*64*32,1200]) ;
    %resb = eegfiltfft(resb,1/.41,0,1) ; 
    %allb = reshape(resb,[2,64,64,32,1200]) ;
    %%% get the triggers and canonical HRF convolved ideal time series
    %boldTR = .275 ; 
    bold_hrf = spm_hrf(boldTR) ; % spm hrf with a tr of 2 seconds
    bold_hrf(2:5) = -.03 ;
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
  
    % *********** bold corrvols *************
    clear boldcorrvols
    boldcorrvol = zeros(1,1,1,size(allb,5)) ; boldcorrvol(1,1,1,:) = boldidealhrf ; 
    boldcorrvol = repmat(boldcorrvol,[size(allb,2),size(allb,3),size(allb,4),1]) ; 
    for i=1:size(allb,1)  
        boldcorrvols(i,:,:,:) = voxcorr(squeeze(allb(i,:,:,:,:)),boldcorrvol) ; 
    end
    postmask = zeros(size(boldcorrvols)) ; postmask(:,1:end,35:end,1:end) = 1 ; 
    boldcorrvols = boldcorrvols.*postmask ; 
    
    bothcorrvols = squeeze((mean(boldcorrvols,1))) ; 

    
    for corr=.08:.02:.45
    
    clear boldepochs 
    for sc=1:size(allb,1) %  for all scans
        inds = find(bothcorrvols>corr & bothcorrvols<corr+.05) ; 
        [i1,i2,i3] = ind2sub(size(bothcorrvols),inds) ; 
        clear boldcorrvoxels boldcorrts
        for i=1:max(size(i1)) ; % get the r values at each index
            boldcorrvoxels(i) = bothcorrvols(i1(i),i2(i),i3(i)) ; 
            boldcorrts(i,:) = allb(sc,i1(i),i2(i),i3(i),:) ; 
        end     
        [cvrvals,cvsortinds] = sort(boldcorrvoxels,'descend') ; 
        sortvox = boldcorrts(cvsortinds,:) ; 
        % get the stimulus indices
        boldtrigvols = allboldtrvols(boldtrigs(sc),:) ; % stimulus onset for that scan
        boldtrigtypes = types(boldtrigs(sc),:) ; 
        for i=1:max(size(stimtypes))
           boldtypinds = boldtrigvols(find(boldtrigtypes==stimtypes(i))) ;  
           for j=1:max(size(boldtypinds))
                boldepochs(i,(sc*10-10)+j,:,:) = sortvox(:,boldtypinds(j)-round(4./boldTR):boldtypinds(j)+round(16./boldTR)) ;                
           end
        end
    end
    
    % baseline correct the single voxel epochs
    clear boldnormepochs
    boldepochs = double(boldepochs) ;  
    for i=1:size(boldepochs,1) 
        for j=1:size(boldepochs,2) ; 
            for k=1:size(boldepochs,3)
                baseijk = squeeze(mean(boldepochs(i,j,k,round(3./boldTR):round(4./boldTR)),4)) ; 
                boldnormepochs(i,j,k,:) = (boldepochs(i,j,k,:)-baseijk)./baseijk ;                
            end            
        end
    end   
    figure,
    
    a = squeeze(mean(mean(boldnormepochs,1),2)) ; 
    zs = zscore(std(a,0,2)) ; goods = find(zs < 0) ; 
    avgvoxtrials = squeeze(mean(boldnormepochs(:,:,goods,:),3)) ;
    allavgvox(sub,:,:,:) = avgvoxtrials ; 
    timess = ((-round(4./boldTR):round(16./boldTR))).*boldTR ;
    for i=1:size(avgvoxtrials,3)
       p1(i) = anova1(squeeze(avgvoxtrials([1,2],:,i))',[],'off') ; 
       p2(i) = anova1(squeeze(avgvoxtrials([1,3],:,i))',[],'off') ; 
    end
    
    %for i=1:3 ; subplot(2,2,i) ; imagesc(squeeze(avgvoxtrials(i,:,:))) ; end
    
    
    ntrials = 20 ; 
    meant = squeeze(mean(avgvoxtrials(:,:,:),2)) ; maxt = max(max(meant)) ; mint = min(min(meant)) ;
    subplot(2,2,1) ; 
    shadedErrorBar([],squeeze(mean(avgvoxtrials(1,:,:),2)),squeeze(std(avgvoxtrials(1,:,:),0,2))./sqrt(ntrials),{'Color',[0,0,1]},0) ; hold on ; ylim([mint-.005,maxt+.01]) ;
    shadedErrorBar([],squeeze(mean(avgvoxtrials(2,:,:),2)),squeeze(std(avgvoxtrials(2,:,:),0,2))./sqrt(ntrials),{'Color',[0,1,0]},0) ; hline(0,'k') ; vline(round(4./boldTR),'r') ; vline(round(6./boldTR),'r') ; 
    xlim([1,size(avgvoxtrials,3)]) ;set(gca,'XTick',1:ceil(size(timess,2)/10):size(timess,2),'XTickLabel',timess(1:ceil(size(timess,2)/10):size(timess,2))) ; xlabel('time(s)') ; title('unperturbed(blue) vs 5% contrast(green)') ; 
    subplot(2,2,2) ; 
    shadedErrorBar([],squeeze(mean(avgvoxtrials(1,:,:),2)),squeeze(std(avgvoxtrials(1,:,:),0,2))./sqrt(ntrials),{'Color',[0,0,1]},0) ; hold on ;  ylim([mint-.005,maxt+.01]) ;
    shadedErrorBar([],squeeze(mean(avgvoxtrials(3,:,:),2)),squeeze(std(avgvoxtrials(3,:,:),0,2))./sqrt(ntrials),{'Color',[0,1,0]},0) ; hline(0,'k') ; vline(round(4./boldTR),'r') ; vline(round(6./boldTR),'r') ; 
    xlim([1,size(avgvoxtrials,3)]) ; set(gca,'XTick',1:ceil(size(timess,2)/10):size(timess,2),'XTickLabel',timess(1:ceil(size(timess,2)/10):size(timess,2))) ; xlabel('time(s)') ; title('unperturbed(blue) vs 60% randomization(green)') ;
    subplot(2,2,3) ; plot(p1,'LineWidth',2,'Color',[0,0,0]) ; xlim([1,size(avgvoxtrials,3)])
    set(gca,'XTick',1:ceil(size(timess,2)/10):size(timess,2),'XTickLabel',timess(1:ceil(size(timess,2)/10):size(timess,2))) ; xlabel('time(s)') ; hline(0.05,'b') ; vline(round(4./boldTR),'r') ; vline(round(6./boldTR),'r') ;
    ylabel('p-value') ; title('p(unperturbed != 5%contrast)') ; 
    subplot(2,2,4) ; plot(p2,'LineWidth',2,'Color',[0,0,0]) ; xlim([1,size(avgvoxtrials,3)])
    set(gca,'XTick',1:ceil(size(timess,2)/10):size(timess,2),'XTickLabel',timess(1:ceil(size(timess,2)/10):size(timess,2))) ; xlabel('time(s)') ; hline(0.05,'b') ; vline(round(4./boldTR),'r') ; vline(round(6./boldTR),'r') ;
    ylabel('p-value') ; title('p(unperturbed != 60%randomization)') ; 
    suptitle(stitle) ; 
    
    end
 