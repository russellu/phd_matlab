cd c:/shared/MONG_01_RB/ ; close all ; clear all ; 
ls 
mongs = dir('MONG_01_RB_RUSSELL_*vhdr') ; 
for i=1:length(mongs) ; 
   EEG = pop_loadbv('.',mongs(i).name) ;  
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
   EEG = denoise_grad(EEG) ; 
   EEG = pop_resample(EEG,512) ; 
   eegs{i} = EEG ; 
   EEG = denoise_bcg(EEG) ;   
   if i==1 ; merged = EEG  ;else merged = pop_mergeset(EEG,merged) ; end     
end
merged = pop_resample(merged,256) ; 
filtmerged = merged ; filtmerged.data = eegfiltfft(filtmerged.data,filtmerged.srate,1,128) ; 
ica = pop_runica(filtmerged,'runica') ; 
for i=1:64 ; subplot(5,13,i) ; topoplot(ica.icawinv(:,i),EEG.chanlocs) ; title(i) ; end
trigs = {'S  1','S  2','S  3','S  4','S  5','S  6'} ;
clear ersp ;
for i=1:size(trigs,2) 
    ep = pop_epoch(ica,{trigs{i}},[-.85,2.85]) ;
    for j=1:64 ; 
        [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(j,:,:)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,'plotersp','off','plotitc','off',...
            'freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',0) ;         
    end
end
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(ersp([1,4],i,:,:),1)),[-4,4]) ; title(i) ; end

% process the FMRI
% get the triggers for each volume separately from the EEG files
cd('c:\shared\MONG_01_RB\mong_rb') ; ls ; clear triginds volOnsets lats
for i=1:length(eegs) ; 
    volOnsets{i} = find(strcmp('R128',{eegs{i}.urevent.type})) ;    
    latsi = cell2mat({eegs{i}.urevent.latency}) ; 
    lats{i} = latsi(volOnsets{i}) ; 
    for t=1:length(trigs)
        trigindst = find(strcmp(trigs{t},{eegs{i}.urevent.type})) ; % trigger indices of stimuli
        triglats = latsi(trigindst) ; % latencies for stimuli
        for lt=1:length(triglats) % for all trials of that stimulus type
           latdiffs = abs(triglats(lt) - latsi(volOnsets{i})) ; % the difference between that trial's latency and all TR latencies
           pnlatdiffs = (triglats(lt) - latsi(volOnsets{i})) ; % pos and negative latdiffs
           triginds(i,t,lt) = find(latdiffs==min(latdiffs)) ; % find the minimum difference
           offsets(i,t,lt) = pnlatdiffs(find(latdiffs==min(latdiffs))) ; % save the offsets from that TR to the actual trigger (positive difference means 
        end
    end
end

% get the FMRI data ; 
TR = 0.9 ; task = round(10/TR) ; 
hrf = spm_hrf(TR) ; 
clear epochs
for i=1:size(triginds,1) ; 
    f1 = load_untouch_nii(['reg_',num2str(i),'.nii.gz']) ; 
    fimg = f1.img ; 
    for j=1:size(triginds,2) ; disp(j) ; 
        for k=1:size(triginds,3)
            epochs(j,(i*10-10)+k,:,:,:,:) = fimg(:,:,:,triginds(i,j,k):triginds(i,j,k)+task) ; 
        end
    end
end


corrs = zeros(size(epochs,1),size(epochs,2),size(epochs,3),size(epochs,4),size(epochs,5)) ; 
hrf = hrf(1:task+1) ; 
for i=1:size(epochs,1) ; disp(i) ; 
    for j=1:size(epochs,2)
        corrs(i,j,:,:,:) = voxcorr(squeeze(epochs(i,j,:,:,:,:)),hrf) ;        
    end
end

f1 = load_untouch_nii('f1_1.nii.gz') ; 
for i=1:6 ; f1.img = ((squeeze(mean(corrs(i,:,:,:,:),2))).^2) ; save_untouch_nii(f1,['corrs_',num2str(i),'.nii.gz']) ; end



mcorrs = squeeze(mean(corrs,2)) ; 
for i=1:36 ; subplot(6,6,i) ;
    
   plotoverlayIntensity2D(squeeze(f1.img(:,:,i)),squeeze(mat2gray(mcorrs(1,:,:,i))),squeeze(mcorrs(1,:,:,i)),270) ;  
    
end


clear rgb ; icount = 1 ; 
for i=7:15
    subplot(3,3,icount) ; 
    img = uint8(mat2gray(double(imrotate(squeeze(f1.img(:,:,i)),270)))*255) ; 
    rgb(:,:,1) = double(imrotate(squeeze(mcorrs(3,:,:,i))>.25,270)) ; 
    rgb(:,:,3) = double(imrotate(squeeze(mcorrs(4,:,:,i))>.25,270)) ; 
    mask=mat2gray(squeeze(sum(rgb,3)==1)) ; 
    imshow(img, 'InitialMag', 'fit') ;
    hold on ; h = imshow(rgb) ; hold off ; 
    set(h, 'AlphaData', (squeeze(mask))) ;
    icount = icount + 1 ; 
end






















