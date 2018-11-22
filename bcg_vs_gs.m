clear all ; close all; 
subs = {'alex','dina','genevieve','jeremie','russell','tegan','valerie'}; 
sets = {'retino_allstim*01*set','retino_allstim*02*set','retino_gamma*01*set','retino_gamma*02*set','retino_movie*set','retino_rest*set'}; 
fmris = {'bp_clean_retino_allstims_01.nii.gz','bp_clean_retino_allstims_02.nii.gz','bp_clean_retino_gamma_01.nii.gz',...
    'bp_clean_retino_gamma_02.nii.gz','bp_clean_retino_movie.nii.gz','bp_clean_retino_rest.nii.gz'};
d1s = {'*retino_allstim*01*1D','*retino_allstim*02*1D','*retino_gamma*01*1D','*retino_gamma*02*1D','*retino_movie*1D','*retino_rest*1D'}; 

fmris = {'bp_gsr_clean_retino_allstims_01.nii.gz','bp_gsr_clean_retino_allstims_02.nii.gz','bp_gsr_clean_retino_gamma_01.nii.gz',...
    'bp_gsr_clean_retino_gamma_02.nii.gz','bp_gsr_clean_retino_movie.nii.gz','bp_gsr_clean_retino_rest.nii.gz'};

freqs = 1:2:100; 

for sb=1:length(subs)
        
    % get the bcg (and electrode bcg) 
    for st=1:length(sets)
       cd(['E:\badger_eeg\',subs{sb}]);
       set_st = dir(sets{st}); 
       eeg = pop_loadset(set_st.name); 
       
       trigs = {eeg.urevent.type};
       r128s = find(strcmpi(trigs,'R128')); 
       lats = cell2mat({eeg.urevent.latency}); 
       r128lats = lats(r128s);     
       % filter the bcg and downsample to number of BOLD samples 
       
       bcgica = load('bcgica'); bcgica = bcgica.bcgica; 
       weights = bcgica{1}; sphere = bcgica{2}; 
       acts = weights*sphere*eeg.data; 
       %eeg.data = acts(1:5,:); 
       eeg.data = eeg.data(32,:); 
       winvs(:,:,sb) = pinv(weights*sphere); 
       filts = zeros(length(freqs),size(eeg.data,1),size(eeg.data,2)); 
       for f=1:length(freqs)
           filts(f,:,:) = abs(eegfiltfft(eeg.data(:,:),eeg.srate,freqs(f)-1,freqs(f)+1)); 
       end
       filts = filts(:,:,r128lats(1):r128lats(end)); 
       clear resfilts convedfilts
       for i=1:length(freqs)
          resfilts(i,:,:) = eegfiltfft(imresize(squeeze(filts(i,:,:)),[size(eeg.data,1),length(r128lats)]),1/0.693,0.01,2);  
       end
       %{
       hrf = spm_hrf(0.693); 
       for i=1:length(freqs)
          for j=1:5
             conved = conv(squeeze(resfilts(i,j,:)),hrf,'full'); 
             conved = conved(1:size(resfilts,3)); 
             convedfilts(i,j,:) = conved; 
          end
       end
        %}
       % get the BOLD
       cd(['E:\sim\fmri\',subs{sb},'\']);
       %bold = load_untouch_nii(fmris{st}); 
       gs_s = dir(d1s{st}); 
       ts = load(gs_s.name); 
       
       ts =  eegfiltfft(ts',1/0.693,0.01,1); 
       
       for i=1:50
           for j=1:size(eeg.data,1)
               xcorrs(sb,st,i,j,:) = xcorr(squeeze(ts(30:size(resfilts,3)-30)),smooth(squeeze(resfilts(i,j,30:end-30))),25,'coeff'); 
           end
       end

    end

end

hrf = spm_hrf(0.693); 
zhrf = zeros(1,51) ; zhrf(26:end) = hrf(1:26); 
bcg_comp = 1; 
times = [-25:25]*0.693; 
subplot(1,3,1); 
imagesc(times,1:2:100,squeeze(mean(xcorrs(:,6,:,bcg_comp,:),1)),[-.15,.15]) ; axis xy ; colormap jet; h = colorbar ; title(h,'corr.'); 
xlabel('time lag (s)') ;ylabel('frequency(hz)'); 
subplot(1,3,2); 
shadedErrorBar(times,squeeze(mean(mean(xcorrs(:,6,4:8,bcg_comp,:),1),3)),squeeze(std(mean(xcorrs(:,6,4:8,bcg_comp,:),3),0,1))/sqrt(7),{'Color',[0,0,1]}); 
hold on ; plot(times,zhrf/2,'k','LineWidth',2); vline(0,'k') ; vline(4.9,'r'); xlabel('time lag (s)'); ylabel('alpha-BOLD correlation'); 

subplot(1,3,3);
plot(1,'b','LineWidth',4); hold on ; plot(1,'k','LineWidth',4);  legend({'BCG alpha xcorr with GS','canonical HRF'});


