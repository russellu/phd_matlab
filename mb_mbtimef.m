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
    size(allb) 
    disp('filtering data') ; 
    %resb = reshape(allb,[2*64*64*32,1200]) ;
    %resb = eegfiltfft(resb,1/.41,0,1) ; 
    %allb = reshape(resb,[2,64,64,32,1200]) ;
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

fb = zeros(size(allb,1),size(allb,2),size(allb,3),size(allb,4),size(allb,5)) ; 
for i=1:size(allb,1)
    for j=1:size(allb,2) ; disp(j) ;
        for k=1:size(allb,3)
            for el=1:size(allb,4)
                fb(i,j,k,el,:) = abs(fft(squeeze(allb(i,j,k,el,:)))) ; 
            end
        end
    end
end
Fs  = 1/.41 ; NFFT = 1200 ;
f = Fs/2*linspace(0,1,NFFT/2+1);
f = f(1:600) ; 
rfb = sqrt(squeeze(mean(fb(:,:,:,:,1:600),1))) ; 
heartinds = find(f<1&f>.8) ; 
heartbs = squeeze(mean(rfb(:,:,:,heartinds),4)) ; 
stiminds = find(f>.055&f<.07) ; 
stimbs = squeeze(mean(rfb(:,:,:,stiminds),4)) ; 
clear rgb 
rgb(:,:,:,1) = heartbs*1.4 ; 
rgb(:,:,:,3) = stimbs ; 
for i=1:32 ; subplot(4,8,i) ; imshow(uint8(mat2gray(squeeze(rgb(:,:,i,:)))*255)) ; end

% do paired sample t test from baseline to task.






%{

%%% get the stimulus spectrogram and plot the power in each voxel.
timess = ((-round(4./boldTR):round(16./boldTR))).*boldTR ;
voxels = [1.17*10^5:1.215*10^5] ; 
resb = reshape(allb,[2*64*64*32,1200]) ;
clear ersp
for i=1:size(resb,1)
    [ersp(i,:,:),itc,powbase,times,freqs,erspboot,itcboot] = newtimef(resb(i,:),size(resb,2),[timess(1),timess(max(size(timess)))],1/.41,0,'plotitc','off','plotersp','off','baseline',NaN,'nfreqs',64) ;
end

ffts = zeros(max(size(resb)),1200) ; 
for i=1:size(resb,1) ;
    disp(num2str(i/262144)) ; 
    ffts(i,:) = abs(fft(squeeze(resb(i,:)))) ;  
end

Fs  = 1/.41 ; NFFT = size(ffts,2) ;
f = Fs/2*linspace(0,1,NFFT/2+1);
f = f(1:600) ; 
logabs = sqrt(ffts(:,1:600)) ; 
incrs = 1:10:600 ; clear meanlogs ; 
for i=1:max(size(incrs))-1
    meanlogs(:,i) = squeeze(mean(logabs(:,incrs(i):incrs(i)+10),2)) ; 
end

for i=1:size(meanlogs,2)
    
    allbs(i,:,:,:) = squeeze(mean(reshape(squeeze(meanlogs(:,i)),[2,64,64,32]),1)) ;
    
end

for i=1:59 ; subplot(6,10,i) ; imagesc(squeeze(allbs(i,:,:,12)),[3,8]) ; end

for fg=1:5:59
    figure,
    
for i=1:32 ; subplot(4,8,i) ; imagesc(sqrt(squeeze(mean(allbs(fg,:,:,i),1))),[10,60]) ; end
end
clear allbs ; 
for i=1:size(logabs,2)  
    allbs(:,:,:,i) = squeeze(mean(reshape(squeeze(logabs(:,i)),[2,64,64,32]),1)) ;   
end

heartinds = find(f<1&f>.8) ; 
heartbs = squeeze(mean(allbs(:,:,:,heartinds),4)) ; 

figure,for i=1:32 ; subplot(6,6,i) ; imagesc(squeeze(heartbs(:,:,i))) ; end
for i=1:32 ; subplot(6,6,i) ; imagesc(squeeze(mean(mean(allb(:,:,:,i,:),1),5))) ; end

save_nii(make_nii(rfb),'rfb.nii.gz') ; 


%}
