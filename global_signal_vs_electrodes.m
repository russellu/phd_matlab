clear all ; close all ; 

subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 

for sb=1:length(subs)
    
    cd(['E:\badger_eeg\',subs{sb}]);
    
    ls('*set'); 
    
    eeg = pop_loadset('retino_rest.set');
    trigs = {eeg.urevent.type};
    trtrigs = find(strcmpi('R128',trigs)); 
    lats = {eeg.urevent.latency}; 
    trlats = cell2mat(lats(trtrigs)); 
    
    %
    freqs = 1:2:100; 
    for elec=1:64
        fcount = 1; 
        filtdat = zeros(50,size(eeg.data,2)); 
        for f=1:length(freqs)
            filtdat(fcount,:) = eegfiltfft(eeg.data(elec,:),eeg.srate,freqs(f)-2,freqs(f)+2);    
            fcount = fcount + 1; 
        end
        filtdat = filtdat(:,trlats(1):trlats(end)); 
        o_logabsfilts = log(abs(filtdat)); 
        o_absfilts = (abs(filtdat)); 
        
        filtdat = imresize(filtdat,[50,length(trlats)]); 
        logabsfilts(elec,:,:) = imresize(o_logabsfilts,[50,length(trlats)]); 
        absfilts(elec,:,:) = imresize(o_absfilts,[50,length(trlats)]); 
        allfilts(elec,:,:) = filtdat; 
    end
    
    for i=1:size(logabsfilts,1)
        for j=1:size(logabsfilts,2)
           logabsij = squeeze(logabsfilts(i,j,:)); 
           goods = find(~isinf(logabsij) & ~isnan(logabsij)); 
           mean_interp = mean(logabsij(goods)); 
           bads = find(isnan(logabsij) | isinf(logabsij)); 
           logabsij(bads) = repmat(mean_interp,[1,length(bads)]); 
           logabsfilts(i,j,:) = logabsij; 
        end
    end
    
    
    res_raweeg = imresize(eeg.data(:,trlats(1):trlats(end)),[64,length(trlats)]); 
    postelecs = [15,51,7,37,19,38,8,52,16,49,45,31,46,60,9,20,10];

    post_raw = squeeze(mean(res_raweeg(postelecs,:),1)); 
    post_absfilts = squeeze(mean(absfilts(postelecs,:,:),1)); 
    post_logabsfilts = squeeze(mean(logabsfilts(postelecs,:,:),1)); 

    all_post_raw(sb,:) = post_raw;
    all_post_absfilts(sb,:,:) = post_absfilts;
    all_post_logabsfilts(sb,:,:) = post_logabsfilts; 
    
    for i=1:50
       abs_xcorr(sb,i,:) = xcorr(post_raw(20:end-20),post_absfilts(i,20:end-20),20,'coeff');          
       logabs_xcorr(sb,i,:) = xcorr(post_raw(20:end-20),post_logabsfilts(i,20:end-20),20,'coeff');  
    end
    
    cd(['E:\fmris\badger_',subs{sb},'\epiregs']);
    fmri = load_untouch_nii('bp_reg_topup_mc_retino_rest.nii.gz'); 
    fmri.img = fmri.img(:,:,:,1:length(trlats)); 
    mask = load_untouch_nii('../f1.nii.gz'); 
    resfmri = reshape(fmri.img,[numel(fmri.img(:,:,:,1)),size(fmri.img,4)]); 
    
    globinds = find(mask.img>50);
    all_global = squeeze(mean(resfmri(globinds,:),1)); 

    for i=1:64
       raw_xcorrs(sb,i,:) = xcorr(all_global(20:end-20),res_raweeg(i,20:end-20),25,'coeff');  
    end
    
    post_rawx(sb,:) = xcorr(all_global(20:end-20),post_raw(20:end-20),20,'coeff'); 
    for i=1:50
       post_absx(sb,i,:) = xcorr(all_global(20:end-20),post_absfilts(i,20:end-20),20,'coeff');  
       post_logabsx(sb,i,:) = xcorr(all_global(20:end-20),post_logabsfilts(i,20:end-20),20,'coeff');       
    end
    
    % bcg correlations
    bcg_rawx(sb,:) = xcorr(all_global(20:end-20),squeeze(res_raweeg(32,20:end-20)),20,'coeff'); 
    for i=1:50
       bcg_absx(sb,i,:) = xcorr(all_global(20:end-20),squeeze(absfilts(32,i,20:end-20)),20,'coeff');  
       bcg_logabsx(sb,i,:) = xcorr(all_global(20:end-20),squeeze(logabsfilts(32,i,20:end-20)),20,'coeff');       
    end
    
    for i=1:64
        for j=1:50
             alles_absx(sb,i,j,:) = xcorr(all_global(20:end-20),squeeze(absfilts(i,j,20:end-20)),20,'coeff');  
             alles_logabsx(sb,i,j,:) = xcorr(all_global(20:end-20),squeeze(logabsfilts(i,j,20:end-20)),20,'coeff');       
        end
    end
    
    
    %{
    brainvox = resfmri(globinds,:); 
    vox_corrs = zeros(64,size(brainvox,1),51); 
    for i=1:64 ; disp(i);
        for j=1:size(brainvox,1)             
            vox_corrs(i,j,:) = xcorr(brainvox(j,20:end-20),res_raweeg(i,20:end-20),25,'coeff'); 
        end
    end
    
    href = load_untouch_nii('../hr_ref.nii.gz'); 
    zref = zeros(size(href.img)); 
    elabs ={eeg.chanlocs.labels};
    %}
    
    %}
    
    %{
    cd(['E:\fmris\badger_',subs{sb}]);
    mkdir xcorr_elecs; 
    cd xcorr_elecs; 
    for elec=1:64
        for i=1:51
            m_voxcorrs = squeeze(mean(vox_corrs(elec,:,i),1)); 
            zmask = zeros(size(mask.img));     
            zmask(globinds) = m_voxcorrs; 
            zref(:,:,:,i) = zmask;
        end
        href.img = zref ; save_untouch_nii(href,['rest_',elabs{elec},'_bcg_vs_bold.nii.gz']); 
    end
    
    
    subplot(2,2,1) ; 
    shadedErrorBar(times,squeeze(mean(raw_xcorrs(:,32,:))),squeeze(std(raw_xcorrs(:,32,:),0,1))/sqrt(8)); 
    ylim([-.2,.1]);hline(0,'k'); vline(0,'k'); title('bcg'); 
    
    subplot(2,2,2) ; 
    shadedErrorBar(times,squeeze(mean(mean(raw_xcorrs(:,[1:31,33:end],:),1),2)),squeeze(std(mean(raw_xcorrs(:,[1:31,33:end],:),2),0,1))/sqrt(8)); 
    ylim([-.2,.1]);hline(0,'k'); vline(0,'k'); title('all channels'); 
    
    subplot(2,2,3); 
    topoplot(squeeze(mean(raw_xcorrs(:,:,25),1)),eeg.chanlocs,'electrodes','labels');

    subplot(2,2,4) ; 
    shadedErrorBar(times,squeeze(mean(mean(raw_xcorrs(:,[46],:),1),2)),squeeze(std(mean(raw_xcorrs(:,[46],:),2),0,1))/sqrt(8)); 
    ylim([-.3,.25]);hline(0,'k'); vline(0,'k'); title('PO4');
    
    
    for i=1:64 
        subplot(5,13,i);
        shadedErrorBar([],squeeze(mean(mean(raw_xcorrs(:,i,:),1),2)),squeeze(std(mean(raw_xcorrs(:,i,:),2),0,1))/sqrt(8)); title(i); 
    end
    topoplot(squeeze(mean(raw_xcorrs(:,:,25),1)),eeg.chanlocs);
    
    %}
        
end
    
%{
times = 0.693*(-25:25); 
subplot(1,2,1)
imagesc(times,1:2:100,squeeze(mean(post_absx,1)),[-.15,.15]) ; axis xy ; colormap jet; vline(0,'k') ; title('abs'); 
subplot(1,2,2); 
imagesc(times,1:2:100,squeeze(mean(post_logabsx,1)),[-.15,.15]) ; axis xy ; colormap jet; vline(0,'k'); title('log'); 

subplot(1,2,1)
imagesc(times,1:2:100,squeeze(mean(abs_xcorr,1)),[-1,1]) ; axis xy ; colormap jet; vline(0,'k') ; title('abs'); 
subplot(1,2,2); 
imagesc(times,1:2:100,squeeze(mean(logabs_xcorr,1)),[-1,1]) ; axis xy ; colormap jet; vline(0,'k'); title('log'); 

%}



%{
allniis = zeros(64,64,33,51,8); 
clear allniis; 
for sb=1:length(subs)
    cd(['E:\fmris\badger_',subs{sb},'\xcorr_elecs']); 
    ls 
    for elec=1%:length(elabs)
       nii = load_untouch_nii(['rest_',elabs{32},'_bcg_vs_bold.nii.gz']); 
       allniis(:,:,:,:,1,sb) = nii.img; 
       disp([subs{sb},' ',elec]); 
    end
end
allniis(isnan(allniis(:))) = 0; 
t_0 = (squeeze(mean(allniis(:,:,:,25,1,:),6)));
t_1 = (squeeze(mean(allniis(:,:,:,32,1,:),6)));

icount=1;for i=6:15 ; subplottight(5,5,icount) ; imagesc(imrotate(squeeze(t_0(:,:,i)),270),[-.15,.15]) ; icount = icount + 1 ; colormap jet; end

icount=1;for i=6:15 ; subplottight(5,5,icount+15) ; imagesc(imrotate(squeeze(t_1(:,:,i)),270),[-.15,.15]) ; icount = icount + 1 ; colormap jet; end


icount=1;for i=1:2:50 ; subplottight(5,5,icount) ; imagesc(imrotate(squeeze(mean(allniis(:,:,12,i,1,:),6)),270),[-.15,.15]) ; colormap jet; icount=icount+1; end

%}

