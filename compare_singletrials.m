clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie'} ; 
eegcomps = {[30,26],[32,23],[21,28],[47,46],[19,37],[47,38]} ;

for sub=1%:length(subs)
    %%%% BOLD processing:
    cd(['c:/shared/badger_mri/',subs{sub},'/nii/melodic']) ; 
    mix = load('melodic_mix') ; 
    weights = load_untouch_nii('melodic_IC.nii.gz') ; 
    allweights{sub} = weights.img ; 
    meansub = load_untouch_nii('mean.nii.gz') ; 
    allmeans{sub} = meansub.img ; 
    mixinds = 1:735:735*5 ; 
    for i=1:length(mixinds) ; boldinds{i} = mixinds(i):mixinds(i)+734 ; end
    boldinds{6} = 735*5+1:735*5+1+449 ; 
    for i=1:length(boldinds) ; segmix{i} = mix(boldinds{i},:) ; end
    bold1 = load_untouch_nii('../reg_topup_mc_retino_gamma_01.nii.gz') ; 
    bold2 = load_untouch_nii('../reg_topup_mc_retino_gamma_02.nii.gz') ;  
    
    boldimgs(1,:,:,:,:) = bold1.img ; boldimgs(2,:,:,:,:) = bold2.img ; 
    
    %%%% EEG processing:
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    prefix = 'highfreq_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ; gammas2=dir([prefix,'preproc*gamma*set']) ;
    setNames = {gammas1(1).name,gammas2(2).name} ;   
    finalersp = zeros(3,64,50,200) ; 
    %figure,
    allmeanepochs = zeros(3,16,16) ; 
    allbersp = zeros(3,50,100) ; 
    for eset=1:2
        EEG = pop_loadset(setNames{eset}) ; 
        eegtrigs = {EEG.urevent.type} ; 
        eeglats = cell2mat({EEG.urevent.latency}) ; 
        r128s = find(strcmp(eegtrigs,'R128')) ;
        rlats = eeglats(r128s) ; 

        startT = rlats(1)./EEG.srate ; 

        currentcomps = eegcomps{sub} ; 
        trigs = {'S  1','S  2','S  3'} ; 
        ep = pop_epoch(EEG,{'S  1','S  2','S  3'},[-2,7]) ; 
        for i=1:length(currentcomps)
            for j=1:48
                [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(currentcomps(i),:,j)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                                        'plotersp','off','plotitc','off','baseline',NaN,'timesout',100,'freqs',[1,100],'nfreqs',50,'winsize',64) ; 

            end
        end
        bersp = ersp - repmat(mean(ersp(:,:,:,times<0),4),[1,1,1,100]) ; 
          
        % display the single trials and 
        a = {ep.epoch.eventtype} ; 
        figure,
        for i=1:48 ; 
            subplot(6,8,i) ; 
            ai = a{i} ; 
            trigind = 0 ; 
            for t=1:length(trigs) ; if ~isempty(find(strcmp(trigs{t},ai))) ; trigind = t ; end ; end
            triginds(i) = trigind ; 
            imagesc(squeeze(mean(ersp(:,i,:,:)))) ; axis xy ;
            title(trigs{trigind}) ; 
        end
        
        st1 = find(triginds==1) ; 
        st2 = find(triginds==2) ; 
        st3 = find(triginds==3) ; 
           
        allbersp(1,:,:) = squeeze(allbersp(1,:,:)) + squeeze(mean(bersp(1,st1,:,:),2))./2 ; 
        allbersp(2,:,:) = squeeze(allbersp(2,:,:)) + squeeze(mean(bersp(1,st2,:,:),2))./2 ; 
        allbersp(3,:,:) = squeeze(allbersp(3,:,:)) + squeeze(mean(bersp(1,st3,:,:),2))./2 ; 

        stimbersp(1,(eset-1)*16+1:(eset-1)*16+1+15,:,:) = squeeze(bersp(1,st1,:,:)) ; 
        stimbersp(2,(eset-1)*16+1:(eset-1)*16+1+15,:,:) = squeeze(bersp(1,st2,:,:)) ; 
        stimbersp(3,(eset-1)*16+1:(eset-1)*16+1+15,:,:) = squeeze(bersp(1,st3,:,:)) ; 
        
        allstimbersp((eset-1)*48+1:(eset-1)*48+1+47,:,:) = squeeze(bersp) ; 
        
        % get the latency for each trigger for each stimulus type
        alltrigs = {EEG.urevent.type} ; 
        lats = cell2mat({EEG.urevent.latency}) ; 
        clear stimlats
        for t=1:length(trigs)
            tr = find(strcmp(trigs{t},alltrigs)) ;
            stimlats(t,:) = lats(tr) ; 
        end
        
        stimlats = stimlats./EEG.srate - startT ; 
        stimtrs = round(stimlats./0.693) ; 
        fmricomps = segmix{2+eset} ; 
        fmricomps = eegfiltfft(fmricomps',1/.693,0.02,1.5) ; fmricomps = fmricomps'  ;

        % first, create the ideal time series, convolve, and correlate to
        % isolate components that were modulated by the stimulus presentation
        ideal = zeros(1,735) ; 
        stimduration = round(5./0.693) ; 
        res_stimtrs = reshape(stimtrs,[1,numel(stimtrs)]) ; 
        for i=1:length(res_stimtrs) ; ideal(res_stimtrs(i):res_stimtrs(i)+stimduration) = 1 ; end
        hrf = spm_hrf(0.693) ; 
        conved = conv(ideal,hrf,'full') ; conved = conved(1:735) ;        
        corrs = corr(conved(25:end-25)',fmricomps(25:end-25,:)) ;    

        % plot the ideal beside the BOLD component time series
        %figure,plot(mat2gray(conved)) ; hold on ; plot(mat2gray(fmricomps(:,(corrs==max(corrs))))+.2,'r') ; title(num2str(corrs(corrs==max(corrs)))) ; vline(res_stimtrs,'k') ; 
        compinds = find(corrs==max(corrs)) ; 
       
        % for the voxels 
        corrbrain = voxcorr(squeeze(boldimgs(eset,:,:,:,25:end-25)),conved(25:end-25)) ; 
        si = find(corrbrain>.2) ; 
        sv = corrbrain(si) ; 
        [sortv,sorti] = sort(sv,'descend') ; 
        goodvoxels = si(sorti(1:200)) ; 
        highcorrs = zeros(size(corrbrain)) ; 
        highcorrs(goodvoxels) = 1 ; 
        
        [gx,gy,gz] = ind2sub(size(corrbrain),goodvoxels) ; 
        for i=1:length(goodvoxels)
            ts(i,:) = squeeze(boldimgs(eset,gx(i),gy(i),gz(i),:)) ; 
        end
        
        clear compepochs voxepochs
        for i=1:size(stimtrs,1)
            for j=1:size(stimtrs,2)
                for k=1:length(compinds)
                    compepochs(k,i,j,:) = fmricomps(stimtrs(i,j):stimtrs(i,j)+15,compinds(k)) ;   
                    voxepochs(i,j,:,:) = ts(:,stimtrs(i,j):stimtrs(i,j)+15) ; 
                end
            end
        end
        bvoxepochs = (voxepochs - repmat(voxepochs(:,:,:,1),[1,1,1,size(voxepochs,4)]))./repmat(voxepochs(:,:,:,1),[1,1,1,size(voxepochs,4)]) ; 
        mvoxepochs = squeeze(mean(bvoxepochs,3)) ; 
        allmvoxepochs(:,(eset-1)*16+1:(eset-1)*16+1+15,:) = mvoxepochs ; 
        meants = squeeze(mean(allmvoxepochs(:,:,10:15),3)) ; 
        %subplot(1,2,eset) ; plot(squeeze(mean(mvoxepochs,2))')
        allmeanepochs = allmeanepochs + mvoxepochs./2 ; 
    end 
    
    subplot(2,3,sub) ; plot(squeeze(mean(allmeanepochs,2))')
    subepochs(sub,:,:,:) = allmeanepochs ; 
    subbersp(sub,:,:,:) = allbersp ; 
    
    %%% check the single trial (spontaneous coupling)
    mstimbersp  = squeeze(mean(stimbersp(:,:,:,times>0 & times<5),4)) ; 
    for i=1:50 ; corrs(i) = corr2(meants(3,:),squeeze(mstimbersp(3,:,i))) ; end
    
    
    
    
    
    
end


%{
trtimes = 0:.693:16 ; 
subplot(1,2,1) ; errorbar(squeeze(mean(mean(subepochs,1),3))',squeeze(std(mean(subepochs,3),0,1))'./sqrt(6)) ; legend({'0%rnd','10%rnd','100%rnd'}) ; xlim([0,17]) ; 
ylabel('BOLD %change') ; set(gca,'XTick',1:2:16,'XTickLabel',round(trtimes(1:2:16))) ; xlabel('time(s)') ; vline([1,8]) ; 
subplot(1,2,2) ; barwitherr(squeeze(std(mean(mean(subepochs(:,:,:,10:14),3),4),0,1)./sqrt(6)),squeeze(mean(mean(mean(subepochs(:,:,:,10:14),1),3),4))) ; set(gca,'XTickLabel',{'0%rnd','10%rnd','100%rnd'}) ; 
xlabel('stim type') ; ylabel('BOLD %change') ; 

subplot(2,2,1) ; errorbar(squeeze(mean(mean(subbersp(:,:,:,times>0&times<5),1),4))',squeeze(std(mean(subbersp(:,:,:,times>0&times<5),4),0,1))'./sqrt(6)) ; 
set(gca,'XTick',1:10:50,'XTickLabel',round(freqs(1:10:50))) ; xlabel('frequency(hz)') ; ylabel('power(db)') ; hline(0,'k') ; legend({'0%rnd','10%rnd','100%rnd'}) ; xlim([0,51]) ; title('grand avg power spectrum') ; 
subplot(2,2,2) ; barwitherr(squeeze(std(mean(mean(subbersp(:,:,freqs>30 & freqs<90,times>0&times<5),3),4),0,1))./sqrt(6),squeeze(mean(mean(mean(subbersp(:,:,freqs>30&freqs<90,times>0&times<5),1),3),4))) ; 
xlabel('stim type') ; set(gca,'XTickLabel',{'0%','10%','100%'}) ; title('gamma(30-90Hz)') ; ylabel('power(db)') ; 
subplot(2,2,3) ; barwitherr(squeeze(std(mean(mean(subbersp(:,:,freqs>8 & freqs<25,times>0&times<5),3),4),0,1))./sqrt(6),squeeze(mean(mean(mean(subbersp(:,:,freqs>8&freqs<25,times>0&times<5),1),3),4))) ; 
xlabel('stim type') ; set(gca,'XTickLabel',{'0%','10%','100%'}) ; title('alpha/beta(8-25Hz)') ; ylabel('power(db)') ; 
%}




