clear all ; close all ; 
cd c:/shared/asl2 ; ls 
subjs = {
    'genevieve_multiband_2' 
    'guillaume_multiband_2' 
    'russell_multiband_2' 
} ;
stimtypes = [1,2,3,4,5,6] ; 
sub = 1 ; 
for sub=2%:max(size(subjs)) ;
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
    
    %%% get the triggers and canonical HRF convolved ideal time series
    disp(['calculating correlations for ',subjs{sub}]) ;
    boldtrigs = [2,4] ;  
    boldTR = 0.41 ; 
    %boldTR = 2 ; 
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
  
   % boldstimvols = zeros(2,6,10,64,64,32,42) ;
   clear boldstimvols ; 
    for i=1:max(size(boldtrigs)) 
        typesi = types(boldtrigs(i),:) ; 
        stimcounts = [1,1,1,1,1,1] ;
        for j=1:max(size(typesi)) ;
            stimind = find(stimtypes==typesi(j)) ; 
            timeind = boldtrvols(j) ; 
            boldstimvols(i,stimind,stimcounts(stimind),:,:,:) = (squeeze(mean(allb(i,:,:,:,round(timeind+2/boldTR):round(timeind+5/boldTR)),5))-squeeze(mean(allb(i,:,:,:,round(timeind-1/boldTR):timeind),5))) ./ squeeze(mean(allb(i,:,:,:,round(timeind-1/boldTR):timeind),5)) ; 
            % 
            %(squeeze(mean(allb(i,:,:,:,timeind+3:timeind+12),5))-squeeze(mean(allb(i,:,:,:,timeind-4:timeind),5))) ./ squeeze(mean(allb(i,:,:,:,timeind-4:timeind),5)) ; 
            stimcounts(stimind) = stimcounts(stimind) + 1 ; 
        end
    end
    mboldstims = squeeze(mean(mean(boldstimvols,1),3)) ; 
    tcount = 1 ; clear boldtrials ; 
    for i=1:2
        for j=1:5
            boldtrials(:,tcount,:,:,:) = boldstimvols(i,:,j,:,:,:) ;
            tcount = tcount + 1 ;
        end
    end
   
    %{
   clear mts ;
    for i=2:6 ; disp(i) ; 
        for j=1:42 ;
            mts(i-1,:,:,:,j) = (squeeze(mean(boldtrials(1,:,:,:,:,j),2))-squeeze(mean(boldtrials(i,:,:,:,:,j),2)))./(squeeze(std(boldtrials(1,:,:,:,:,j),0,2))+squeeze(std(boldtrials(i,:,:,:,:,j),0,2))) ;
        end
    end
    %}
    
    
   % boldtrials = squeeze(mean(boldtrials(:,:,:,:,:,round(6/boldTR):round(10/boldTR)),6))-squeeze(mean(boldtrials(:,:,:,:,:,round(2/boldTR):round(3/boldTR)),6)) ;
    
    
    
    
    
    
    mb = squeeze(mean(mean(allb,1),5)) ; 
    sumx = squeeze(sum(sum(mb,2),3)) ; 
    sumy = squeeze(sum(sum(mb,1),3)) ; 
    sumz = squeeze(sum(sum(mb,1),2)) ; 
    
    boldtrials(isnan(boldtrials)) = 0 ; boldtrials(isinf(boldtrials)) = 0 ;
    fvals = zeros(5,size(boldtrials,3),size(boldtrials,4),size(boldtrials,5)) ; 
    for stype = 2:6
    for i=1:size(boldtrials,3) ; 
        disp(i) ; 
        for j=1:size(boldtrials,4) ; 
            for k=1:size(boldtrials,5)
                if sumz(k)>1000
                 [p,anovatab,stats]= anova1(squeeze(boldtrials([1,stype],:,i,j,k))',[],'off') ; 
                 %[h,p,ci,stats ] = ttest(squeeze(boldtrials(1,:,i,j,k)),squeeze(boldtrials(stype,:,i,j,k))) ; 
                 fvals(stype-1,i,j,k) = anovatab{2,5} ; 
                end
            end
        end
    end
    end
    
    clear mts ;
    for i=2:6 ;
        mts(i-1,:,:,:) = (squeeze(mean(boldtrials(1,:,:,:,:),2))-squeeze(mean(boldtrials(i,:,:,:,:),2)))./(squeeze(std(boldtrials(1,:,:,:,:),0,2))+squeeze(std(boldtrials(i,:,:,:,:),0,2))) ;
    end
    
    refimg = load_untouch_nii('slice_mc_singleband_1.nii.gz');
    refimg.img = medfilt3(squeeze(mean(abs(fvals(:,:,:,:)),1))) ; save_untouch_nii(refimg,'fvals.nii.gz') ; 
    mf = squeeze(mean(abs(fvals),1)) ; 
    mmf = medfilt3(mf) ; 
    disp3d(mmf>1.2) ; 
    
    
    %for i=1:10 ; subplot(3,4,i) ; imagesc(squeeze(mean(abs(fvals(:,:,:,6+i)),1))) ; end

   
end





