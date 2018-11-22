clear all ; close all ; 

subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 

for sb=1:length(subs)
    
    cd(['E:\badger_eeg\',subs{sb}]);   
    ls('*set'); 
    
    eeg = pop_loadset('retino_rest.set');
    
    %filtdat = eegfiltfft(eeg.data,eeg.srate,1,120);
    
    %[weights,sphere] = runica(filtdat(:,1:2:end),'maxsteps',128); 
    
    %eeg.data = weights*sphere*eeg.data; 
    postelecs = [15,51,7,37,19,38,8,52,16,49,45,31,46,60,9,20,10];

    trigs = {eeg.urevent.type};
    trtrigs = find(strcmpi('R128',trigs)); 
    lats = {eeg.urevent.latency}; 
    trlats = cell2mat(lats(trtrigs)); 
    
    freqs = 1:2:100;
    elecs_abs = zeros(50,64,length(trlats));
    elecs_logabs = zeros(50,64,length(trlats)); 
    for hz=1:length(freqs)
       filtdata = eegfiltfft(eeg.data,eeg.srate,freqs(hz)-3,freqs(hz)+3);  
       filtdata = filtdata(:,round(trlats(1)):round(trlats(end)));  
       absfiltdata = abs(filtdata); 
       logabsfiltdata = log(abs(filtdata)); 
       
       elecs_abs(hz,:,:) = imresize(absfiltdata,[64,length(trlats)]); 
       elecs_logabs(hz,:,:) = imresize(logabsfiltdata,[64,length(trlats)]); 
        
    end
    
    for i=1:size(elecs_logabs,1)
        for j=1:size(elecs_logabs,2)
           logabsij = squeeze(elecs_logabs(i,j,:)); 
           goods = find(~isinf(logabsij) & ~isnan(logabsij)); 
           mean_interp = mean(logabsij(goods)); 
           bads = find(isnan(logabsij) | isinf(logabsij)); 
           logabsij(bads) = repmat(mean_interp,[1,length(bads)]); 
           elecs_logabs(i,j,:) = logabsij; 
        end
    end
    
    abs_gfp = squeeze(mean(mean(elecs_abs)));
    logabs_gfp = squeeze(mean(mean(elecs_logabs))); 
    
    for i=1:50
       elecs_logabs(i,:,:) = eegfiltfft(squeeze(elecs_logabs(i,:,:)),1/0.693,0.005,2);  
       elecs_abs(i,:,:) = eegfiltfft(squeeze(elecs_abs(i,:,:)),1/0.693,0.005,2);  
    end
     
    cd(['E:\fmris\badger_',subs{sb},'\epiregs']);
    fmri = load_untouch_nii('bp_reg_topup_mc_retino_rest.nii.gz'); 
    fmri.img = fmri.img(:,:,:,1:length(trlats)); 
    mask = load_untouch_nii('c:/shared/epireg/meanmask.nii.gz'); 
    meannii = load_untouch_nii('c:/shared/epireg/mean.nii.gz'); 
    resfmri = reshape(fmri.img,[numel(fmri.img(:,:,:,1)),size(fmri.img,4)]); 
    
    globinds = find(mask.img==1);
    all_global = squeeze(mean(resfmri(globinds,:),1)); 
    
    for i=1:50
        for j=1:64
            logabs_xcorrs(sb,i,j,:) = xcorr(all_global(20:end-20),squeeze(elecs_logabs(i,j,20:end-20)),20,'coeff'); 
            abs_xcorrs(sb,i,j,:) = xcorr(all_global(20:end-20),squeeze(elecs_abs(i,j,20:end-20)),20,'coeff'); 
        end
    end
    
    absgfp_xcorr(sb,:) = xcorr(all_global,abs_gfp,20,'coeff'); 
    logabsgfp_xcorr(sb,:) = xcorr(all_global,logabs_gfp,20,'coeff'); 
    
    for i=1:50
       absgfp_fxcorr(sb,i,:) = xcorr(all_global(20:end-20),squeeze(mean(elecs_abs(i,[1:31,33:end],20:end-20),2)),20,'coeff');  
       logabsgfp_fxcorr(sb,i,:) = xcorr(all_global(20:end-20),squeeze(mean(elecs_logabs(i,[1:31,33:end],20:end-20),2)),20,'coeff');  
    end
    
    hr_xcorrs = zeros(size(globinds,1),41); 
    mean_alpha_hr = squeeze(mean(mean(elecs_logabs(4:7,postelecs,:),1))); 
    conv_alpha_hr = conv(mean_alpha_hr,spm_hrf(0.693),'full');
    conv_alpha_hr = conv_alpha_hr(1:length(mean_alpha_hr)); 
    for i=1:length(globinds)
        hr_xcorrs(i,:) = xcorr(resfmri(globinds(i),20:end-20),mean_alpha_hr(20:end-20),20,'coeff');    
    end
    times = 0.693*(-20:20); 


    conv_hr_xcorr(:,:,:,sb) = voxcorr(fmri.img(:,:,:,20:end-20),conv_alpha_hr(20:end-20)); 
    
    
    for i=1:41
        zimg = zeros(size(mask.img)); 
        zimg(globinds)= hr_xcorrs(:,i); 
      %  subplottight(4,11,i); 
      %  imagesc(imrotate(squeeze(mean(zimg(:,:,8:15),3)),270),[-.2,.2]); title(i); 
      %  plotoverlayIntensity2D(meannii.img(:,:,11),mat2gray(abs(squeeze(mean(zimg(:,:,8:15),3)))),(squeeze(mean(zimg(:,:,8:15),3))),270);
      %  title(['t=',num2str(times(i))]); 
        allzimgs(:,:,:,i,sb) = zimg; 

    end
    
    
    
end
times = 0.693*(-20:20); 


allzimgs(1,1,:,:,:) = -.15; allzimgs(1,2,:,:,:) = .15; 
allzimgs(isnan(allzimgs)) = 0; 
for i=1:41
    subplottight(4,11,i); 
    img = squeeze(mean(mean(allzimgs(:,:,8:15,i,:),3),5)); 
    plotoverlayIntensity2D(meannii.img(:,:,11),mat2gray(abs(img)),img,270);
    title(['t=',num2str(times(i))]); 
end

subplot(1,2,1); 
imagesc(times,1:2:100,squeeze(mean(mean(logabs_xcorrs(:,:,32,:),1),3)),[-.15,.15]) ; colormap jet; axis xy  ; title('bcg channel')
subplot(1,2,2); 
imagesc(times,1:2:100,squeeze(mean(mean(logabs_xcorrs(:,:,postelecs,:),1),3)),[-.15,.15]) ; colormap jet; axis xy  ; xlabel('time lag(s)') ; ylabel('frequency(hz)'); title('posterior channels'); 

%{
imagesc(times,1:2:100,squeeze(mean(mean(abs_xcorrs,1),3)),[-.15,.15]); axis xy ; vline(0,'k'); 
hrf = spm_hrf(0.693); 
hrfvec = zeros(1,41) ; hrfvec(21:end) = hrf(1:21)*0.7; 
plot(squeeze(mean(mean(abs_xcorrs(:,1:4,32,:),1),2)),'LineWidth',2) ; hold on ; plot(squeeze(mean(mean(abs_xcorrs(:,3:5,1:10,:),1),2))','r'); vline(21,'k');
plot(hrfvec,'k','LineWidth',2); 
vline(21,'k'); vline(29,'r'); 
plot(mean(logabsgfp_xcorr),'m','LineWidth',2)
plot(squeeze(mean(logabsgfp_fxcorr(:,6,:),1)),'k'); 
plot(squeeze(mean(mean(logabs_xcorrs(:,5:7,32,:),1),2)),'r','LineWidth',3) ; 
%}

    