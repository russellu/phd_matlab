clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','valerie'} ;
comps = {[32,6,48],[9,18,20],[12,10,4],[24,9,16],[8,15,52],[11,21,33],[46,33,15],[17,7,21]} ;
occ = {[65,13],[44,41,28],[60,48],[85,82,67,6],[81,54,36],[77,76,74],[53,49],[77,75,42]} ; 
lingual = {[80,16],[8],[18],[50,5],[60],[53],[31],[16]} ; 
lateral = {[17],[23],[41],[21],[19],[51],[33],[22]} ; 
threshs = [.25,.25,.25,.25,.25,.25,.25,.25] ;
freqs = 1:100 ;  
clear tcorrs ; 
for sub=1:length(subs)
    cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls 
    bcgicas = dir('bcgica*gamma*set') ; 
    
    clear freqepochs boldepochs
    for scan=1:2 
        cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls 
        EEG = pop_loadset(bcgicas(scan).name) ; 
        EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
        allica = load('highfreqs') ; allica = allica.highfreqs ; 
        weights = allica{1} ; sphere = allica{2} ; 
        ica = EEG ; ica.data = weights*sphere*EEG.data ; 
        winv = pinv(weights*sphere) ;
        compacts = ica.data(comps{sub},:) ; 
        fcomps = zeros(size(compacts,1),length(freqs),size(compacts,2)) ; 
        for i=1:length(freqs)
            fcomps(:,i,:) = eegfiltfft(compacts,250,freqs(i)-4,freqs(i)+4) ; 
        end
        mcomps = squeeze(mean(abs(fcomps))) ; 
        events = {ica.urevent.type} ;
        lats = {ica.urevent.latency} ; lats = cell2mat(lats) ; 
        r128s = find(strcmpi('R128',events)) ; 
        startT = round(lats(r128s(1))) ; endT = round(lats(r128s(end))) ;
        stimlats(1,:) = lats((strcmpi('S  1',events)))-startT ; 
        stimlats(2,:) = lats((strcmpi('S  2',events)))-startT ; 
        stimlats(3,:) = lats((strcmpi('S  3',events)))-startT ; 
        scancomps = mcomps(:,startT:endT) ; 
        ntrs = length(r128s) ; 
        rescomps = imresize(scancomps,[size(scancomps,1),ntrs]) ; 
        res_factor = size(scancomps,2)/ntrs ; 
        stimlats = round(stimlats/res_factor) ; 
        
        %fmri
        cd(['c:/shared/newbadger_mri/',subs{sub}]) ;
        corrs = load_untouch_nii('corrs.nii.gz') ; gamma1 = load_untouch_nii(['warp_gamma_',num2str(scan),'.nii.gz']) ; 
        bincorrs = corrs.img>0.2 ; inds = find(bincorrs==1) ; [cx,cy,cz] = ind2sub(size(corrs.img),inds) ; clear mts ; 
        for i=1:length(inds) ; mts(i,:) = squeeze(gamma1.img(cx(i),cy(i),cz(i),:)) ; end ; mts = mean(mts) ; 
        hrf = spm_hrf(0.693) ; 
        % the resconvs
        padconvs = zeros(size(rescomps,1),size(rescomps,2)+100) ; 
        padconvs(:,101:end) = rescomps ; padconvs(:,1:100) = repmat(mean(rescomps,2),[1,100]) ; 
        for i=1:size(rescomps,1)
            convi = conv(log(padconvs(i,:)),hrf,'full') ;
            padconvs(i,:) = convi(1:size(padconvs,2)) ; 
        end
        padconvs = padconvs(:,101:end) ;
        
        mts(1:10) = 0 ;

        mts = mts(1:size(padconvs,2)) ; 
        for i=1:100
             freqcorrs(sub,scan,i) = corr2(padconvs(i,10:end-10),mts(10:end-10)) ;  
        end
        
        resimg = reshape(gamma1.img,[numel(gamma1.img(:,:,:,1)),735]) ; 
        voxstims = zeros(size(stimlats,1),size(stimlats,2),size(resimg,1),21) ; 

        for i=1:size(stimlats,1)
            for j=1:size(stimlats,2)
                freqepochs(scan,i,j,:,:) = padconvs(:,stimlats(i,j):stimlats(i,j)+20) ; 
                substims(sub,scan,i,j,:,:) = rescomps(:,stimlats(i,j)-5:stimlats(i,j)+12) ; 
                voxstims(i,j,:,:) = resimg(:,stimlats(i,j):stimlats(i,j)+20) ; 
            end
        end
        allvoxstims(scan,:,:,:,:) = voxstims ; 
    end
    allfreqepochs(:,1:16,:,:) = squeeze(freqepochs(1,:,:,:,:)) ; allfreqepochs(:,17:32,:,:) = squeeze(freqepochs(2,:,:,:,:)) ; 
    nvoxstims(:,1:16,:,:) = squeeze(allvoxstims(1,:,:,:,:)) ; nvoxstims(:,17:32,:,:) = squeeze(allvoxstims(2,:,:,:,:)) ; 
    basevox = nvoxstims - repmat(nvoxstims(:,:,:,1),[1,1,1,21]) ; 
    allfreqepochs = allfreqepochs - repmat(mean(allfreqepochs(:,:,:,1),4),[1,1,1,21]) ; 
    
    tvox = squeeze(mean(basevox(:,:,:,12:16),4)) ; 
    tfreq = squeeze(mean(allfreqepochs(:,:,:,12:16),4)) ; 
    for s=1:3
        for i=1:100
            [c,p] = corr(squeeze(tfreq(s,:,i))',squeeze(tvox(s,:,:))) ; 
            cbrains(s,i,:) = c ; 
        end
    end
    
    cat = load_untouch_nii('cat.nii.gz') ;
    for i=1:3 ; 
        zimg = zeros(size(cat.img)) ; 
        for j=1:100
            zimg(:,:,:,j) = reshape(squeeze(cbrains(i,j,:)),size(squeeze(zimg(:,:,:,1)))) ; 
        end
        cat.img = zimg ; 
        save_untouch_nii(cat,['singletrial_',num2str(i),'.nii.gz']) ; 
    end
    
end








