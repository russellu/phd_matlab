clear all ; close all ; 
cd c:/shared/asl2 ; ls 
subjs = {
    'genevieve_multiband_2' 
    'guillaume_multiband_2' 
    'russell_multiband_2' 
} ;
stimtypes = [1,2,3,4,5,6] ; 
sub = 1 ; 
for sub=1:max(size(subjs)) ;
    cd(['c:/shared/multiband/',subjs{sub}]) ;
    %asl=dir('reg_asl*') ; 
    bold=dir('reg*multi*gz') ; 
    
    % get the BOLD
    disp(['loading bold images for ',subjs{sub}]) ;
    clear allb
    for i=1:max(size(bold))
        boldnii = load_untouch_nii(bold(i).name) ;  
        boldimg = boldnii.img ; 
        allb(i,:,:,:,:) = boldimg ; 
    end
    size(allb) 
    disp('filtering data') ; 
    resb = reshape(allb,[2*64*64*32,1200]) ;
    resb = eegfiltfft(resb,1/.41,0,1) ; 
    allb = reshape(resb,[2,64,64,32,1200]) ;
    %%% get the triggers and canonical HRF convolved ideal time series
    disp(['calculating correlations for ',subjs{sub}]) ;
    boldtrigs = [1,3] ;  
    boldTR = 0.41 ; 
   % boldTR = 2 ; 
    bold_hrf = spm_hrf(boldTR) ; % spm hrf with a tr of 2 seconds
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

    refimg = load_untouch_nii('slice_mc_singleband_1.nii.gz');
    refimg.img = single(bothcorrvols) ; save_untouch_nii(refimg,'corrmask.nii.gz') ; 
    
    
    clear boldepochs 
    for sc=1:size(allb,1) %  for all scans
        inds = find(bothcorrvols>.2) ; 
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
                boldepochs(i,(sc*5-5)+j,:,:) = sortvox(:,boldtypinds(j)-round(4./boldTR):boldtypinds(j)+round(16./boldTR)) ;                
           end
        end
    end
    
    
        
    % baseline correct the single voxel epochs
    clear boldnormepochs
    boldepochs = double(boldepochs) ;  
    for i=1:size(boldepochs,1) 
        for j=1:size(boldepochs,2) ; 
            for k=1:size(boldepochs,3)
                baseijk = squeeze(mean(boldepochs(i,j,k,1:round(4./boldTR)),4)) ; 
                boldnormepochs(i,j,k,:) = (boldepochs(i,j,k,:)-baseijk) ;%./baseijk ;                
            end            
        end
    end   
    a = squeeze(mean(mean(boldnormepochs,1),2)) ; 
    zs = zscore(std(a,0,2)) ; goods = find(zs < 1) ; 
    avgvoxtrials = squeeze(mean(boldnormepochs(:,:,goods,:),3)) ;
    allavgvox(sub,:,:,:) = avgvoxtrials ; 
end

timess = ((-round(4./boldTR):round(16./boldTR))).*boldTR ;
tpoints = [find(timess==0),find(abs(timess-2)==min(abs(timess-2)))] ;
for mv=1:3 
    gcf = figure,
    mavg = squeeze(mean(allavgvox(mv,:,:,:),1)) ; 
    subplot(2,2,1) ;
    %errorbar(squeeze(mean(mavg([1,2],:,:),2))',squeeze(std(mavg([1,2],:,:),0,2))'./sqrt(10),'LineWidth',2) ; 
    shadedErrorBar([],squeeze(mean(mavg([1],:,:),2)),squeeze(std(mavg([1],:,:),0,2))./sqrt(10),{'Color',[0,0,1]},0) ;  hold on ; 
    shadedErrorBar([],squeeze(mean(mavg([2],:,:),2)),squeeze(std(mavg([2],:,:),0,2))./sqrt(10),{'Color',[0,1,0]},0) ;  
    xlim([1,size(mavg,3)])
    set(gca,'XTick',1:ceil(size(timess,2)/10):size(timess,2),'XTickLabel',timess(1:ceil(size(timess,2)/10):size(timess,2))) ; vline(tpoints(1),'r') ; vline(tpoints(2),'r') ; xlabel('time(s)') ; ylabel('(task-rest)') ; hline(0,'k') ; 
    title('100%contrast (blue) vs 5%contrast (green)') ; 
    subplot(2,2,2) ; 
    shadedErrorBar([],squeeze(mean(mavg([1],:,:),2)),squeeze(std(mavg([1],:,:),0,2))./sqrt(10),{'Color',[0,0,1]},0) ; hold on ; 
    shadedErrorBar([],squeeze(mean(mavg([6],:,:),2)),squeeze(std(mavg([6],:,:),0,2))./sqrt(10),{'Color',[0,1,0]},0) ;  
    xlim([1,size(mavg,3)])
    set(gca,'XTick',1:ceil(size(timess,2)/10):size(timess,2),'XTickLabel',timess(1:ceil(size(timess,2)/10):size(timess,2))); vline(tpoints(1),'r') ; vline(tpoints(2),'r') ; xlabel('time(s)') ; ylabel('(task-rest)') ; hline(0,'k') ; 
    title('0%rnd (blue) vs 60%rnd (green)'); 
    for i=1:size(mavg,3)
        cps(i) = anova1(squeeze(mavg([1,2],:,i))',[],'off') ;
        rps(i) = anova1(squeeze(mavg([1,6],:,i))',[],'off') ;
    end
    subplot(2,2,3) ; bar(cps) ; hline(0.05,'r') ; title('p(100%contrast - 5%contrast != 0)') ;ylabel('p value') ;  %xlim([1,size(mavg,3)])
    set(gca,'XTick',1:ceil(size(timess,2)/10):size(timess,2),'XTickLabel',timess(1:ceil(size(timess,2)/10):size(timess,2))); vline(tpoints(1),'r') ; vline(tpoints(2),'r') ; xlabel('time(s)')  
    subplot(2,2,4) ; bar(rps) ; hline(0.05,'r') ; title('p(0%random - 60%random != 0)') ;ylabel('p value') ;  %xlim([1,size(mavg,3)])
    set(gca,'XTick',1:ceil(size(timess,2)/10):size(timess,2),'XTickLabel',timess(1:ceil(size(timess,2)/10):size(timess,2))) ; vline(tpoints(1),'r') ; vline(tpoints(2),'r') ; xlabel('time(s)')  
    suptitle(['subject ',num2str(mv)]) ; 
    set( gcf, 'Units', 'normalized', 'Position', [0 0 .35 .45] );
end

allavg = squeeze(mean(allavgvox,3)) ; 


gcf = figure ;
subplot(2,2,1) ; 
shadedErrorBar([],squeeze(mean(allavg(:,[1],:),1)),squeeze(std(allavg(:,[1],:),0,1))./sqrt(3),{'Color',[0,0,1]},0) ;  hold on ; 
shadedErrorBar([],squeeze(mean(allavg(:,[2],:),1)),squeeze(std(allavg(:,[2],:),0,1))./sqrt(3),{'Color',[0,1,0]},0) ;  
xlim([1,size(mavg,3)])
set(gca,'XTick',1:ceil(size(timess,2)/10):size(timess,2),'XTickLabel',timess(1:ceil(size(timess,2)/10):size(timess,2))) ; vline(tpoints(1),'r') ; vline(tpoints(2),'r') ; xlabel('time(s)') ; ylabel('(task-rest)') ; hline(0,'k') ; 
title('100%contrast (blue) vs 5%contrast (green)') ; 
subplot(2,2,2) ;
shadedErrorBar([],squeeze(mean(allavg(:,[1],:),1)),squeeze(std(allavg(:,[1],:),0,1))./sqrt(3),{'Color',[0,0,1]},0) ;  hold on ; 
shadedErrorBar([],squeeze(mean(allavg(:,[6],:),1)),squeeze(std(allavg(:,[6],:),0,1))./sqrt(3),{'Color',[0,1,0]},0) ; 
xlim([1,size(mavg,3)])
set(gca,'XTick',1:ceil(size(timess,2)/10):size(timess,2),'XTickLabel',timess(1:ceil(size(timess,2)/10):size(timess,2))) ; vline(tpoints(1),'r') ; vline(tpoints(2),'r') ; xlabel('time(s)') ; ylabel('(task-rest)') ; hline(0,'k') ; 
title('0%rnd (blue) vs 60%rnd (green)'); 
for i=1:size(mavg,3)
    cps(i) = anova1(squeeze(mavg(:,[1,2],i))',[],'off') ;
    rps(i) = anova1(squeeze(mavg(:,[1,6],i))',[],'off') ;
end
subplot(2,2,3) ; bar(cps) ; hline(0.05,'r') ; title('p(100%contrast - 5%contrast != 0)') ;ylabel('p value') ;  %xlim([1,size(mavg,3)])
set(gca,'XTick',1:ceil(size(timess,2)/10):size(timess,2),'XTickLabel',timess(1:ceil(size(timess,2)/10):size(timess,2))) ; vline(tpoints(1),'r') ; vline(tpoints(2),'r') ; xlabel('time(s)')  
subplot(2,2,4) ; bar(rps) ; hline(0.05,'r') ; title('p(0%random - 60%random != 0)') ;ylabel('p value') ;  %xlim([1,size(mavg,3)])
set(gca,'XTick',1:ceil(size(timess,2)/10):size(timess,2),'XTickLabel',timess(1:ceil(size(timess,2)/10):size(timess,2))) ; vline(tpoints(1),'r') ; vline(tpoints(2),'r') ; xlabel('time(s)')  
suptitle('GRAND AVERAGE') ; 
set( gcf, 'Units', 'normalized', 'Position', [0 0 .35 .45] );










