clear all ; close all ; 
cd c:/shared/asl2 ; ls 
subjs = {
    'genevieve_multiband_2' 
    'guillaume_multiband_2' 
    'russell_multiband_2' 
} ;
stimtypes = [1,2,3,4,5,6] ; 
sub = 1 ; 
for sub=1%:max(size(subjs)) ;
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
    %resb = reshape(allb,[2*64*64*32,1200]) ;
    %resb = eegfiltfft(resb,1/.41,0,1) ; 
    %allb = reshape(resb,[2,64,64,32,1200]) ;
    %%% get the triggers and canonical HRF convolved ideal time series
    disp(['calculating correlations for ',subjs{sub}]) ;
    boldtrigs = [2,4] ;  
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
    
    % get the stimulus epochs
    epochtime = -round(3/boldTR):round(15/boldTR) ; 
    voxepochs = zeros(6,10,size(allb,2),size(allb,3),size(allb,4),max(size(epochtime))) ;
    stimtypecounts = [1,1,1,1,1,1] ;
    for scan=1:max(size(boldtrigs)) ;
        for present=1:size(types(scan,:),2) ;
            %disp(types(scan,present)) ; 
            ctype = types(scan,present) ;
            ctime = boldtrvols(present) ;  
            disp(['ctype = ',num2str(ctype),' ctime = ',num2str(ctime)]) ; 
            voxepochs(scan,stimtypecounts(ctype),:,:,:,:) = squeeze(allb(scan,:,:,:,ctime-round(3/boldTR):ctime+round(15/boldTR))) ; 
            stimtypecounts(ctype) = stimtypecounts(ctype) + 1 ; 
        end
    end
end
