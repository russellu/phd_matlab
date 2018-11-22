clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','valerie'} ;
comps = {[32,6,48],[9,18,20],[12,10,4],[24,9,16],[8,15,52],[11,21,33],[46,33,15],[17,7,21]} ;
occ = {[65,13],[44,41,28],[60,48],[85,82,67,6],[81,54,36],[77,76,74],[53,49],[77,75,42]} ; 
lingual = {[80,16],[8],[18],[50,5],[60],[53],[31],[16]} ; 
lateral = {[17],[23],[41],[21],[19],[51],[33],[22]} ; 
threshs = [.25,.25,.25,.25,.25,.25,.25,.25] ;

for sub=1:length(subs)
    cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls 
    bcgicas = dir('bcgica*gamma*set') ; 
    
    clear mallspecs 
    for scan=1:2 
        EEG = pop_loadset(bcgicas(scan).name) ; 
        EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
        allica = load('highfreqs') ; allica = allica.highfreqs ; 
        weights = allica{1} ; sphere = allica{2} ; 
        ica = EEG ; ica.data = weights*sphere*EEG.data ; 
        winv = pinv(weights*sphere) ;
        events = {ica.urevent.type} ;
        lats = {ica.urevent.latency} ; lats = cell2mat(lats) ; 
        r128s = find(strcmpi('R128',events)); 
        startlat = lats(r128s(1)) ; 
        stimlats(1,:) = lats((strcmpi('S  1',events))) ; 
        stimlats(2,:) = lats((strcmpi('S  2',events))) ; 
        stimlats(3,:) = lats((strcmpi('S  3',events))) ; 
        stimtrs(scan,:,:) = round((stimlats - startlat)/(250*0.693)) ; 
    
        stimtrigs = {'S  1','S  2','S  3'} ; clear ersp alleps
        for s=1:length(stimtrigs)
            ep = pop_epoch(ica,{stimtrigs{s}},[-1,6]) ; 
            alleps(s,:,:,:) = ep.data ; 
        end
        compeps = alleps(:,comps{sub},:,:) ; 

        clear allspecs ; wsize = 125 ; wincr = 5 ; 
        for i=1:size(compeps,1)
            for j=1:size(compeps,2)
                for k=1:size(compeps,4)
                    clear specs ; 
                    M = wsize ; Wc = .1 ; kk = (1:M-1); s = sin(Wc*kk)./ kk ; c0 = [Wc,s]; A = toeplitz(c0);
                    [V,evals] = eig(A); % Only need the principal eigenvector
                    [emax,imax] = max(abs(diag(evals)));
                    w = V(:,end-4:end) ;clear g1 
                    g1(:,1:size(w,2)) = w ; 
                    halfl = round(M/2) ; 
                    chan = squeeze(compeps(i,j,:,k)) ; icount = 1 ; 
                    for ind=halfl+1:wincr:length(chan)-M 
                        datai = chan(ind-halfl:ind+M-halfl-1) ;
                        windowed = repmat(datai,[1,size(w,2)]).*w ; 
                        f = abs(fft(windowed,[],1)) ; 
                        spec = f(1:halfl,:) ; 
                        specs(:,icount,:) = spec ; 
                        icount = icount + 1 ; 
                    end
                    allspecs(i,j,k,:,:,:) = specs ; 
                end
            end
        end
        mallspecs(scan,:,:,:,:) = squeeze(mean(mean(allspecs(:,:,:,:,:,:),2),6)) ; 
    end
    allspectrials(:,1:16,:,:) = squeeze(mallspecs(1,:,:,:,:)) ; allspectrials(:,17:32,:,:) = squeeze(mallspecs(2,:,:,:,:)) ;
    base_spectrials = log(allspectrials) - repmat(mean(log(allspectrials(:,:,:,1:50)),4),[1,1,1,size(allspectrials,4)]) ;
    allsubspecs(sub,:,:,:,:) = allspectrials ;
    allsubbasespecs(sub,:,:,:,:) = base_spectrials ; 

    
    cd(['c:/shared/newbadger_mri/',subs{sub},'/melodic']) ;
    mix = load('melodic_mix') ; cd ..
    mix = eegfiltfft(mix',1/0.693,0.02,1.5)' ; clear gammabold
    gammabold(1,:,:) = mean(mix(735*2+1:735*3,occ{sub}),2) ; 
    gammabold(2,:,:) = mean(mix(735*3+1:735*4,occ{sub}),2) ;
    
    clear allboldepochs 
    for scan=1:2
        corrs = load_untouch_nii('corrs.nii.gz') ; gamma1 = load_untouch_nii(['bp_reg_topup_mc_retino_gamma_0',num2str(scan),'.nii.gz']) ; 
        bincorrs = corrs.img>0.2 ; inds = find(bincorrs==1) ; [cx,cy,cz] = ind2sub(size(corrs.img),inds) ; clear mts ; 
        for i=1:length(inds) ; mts(i,:) = squeeze(gamma1.img(cx(i),cy(i),cz(i),:)) ; end ; mts = mean(mts) ; 
        clear boldepochs
        for i=1:size(stimtrs,2)
            for j=1:size(stimtrs,3)
                boldepochs(i,j,:) = mts(stimtrs(scan,i,j):stimtrs(scan,i,j)+20) ; 
            end
        end
        %boldepochs = boldepochs - repmat(mean(boldepochs(:,:,1:5),3),[1,1,size(boldepochs,3)]) ; 
        allboldepochs(scan,:,:,:) = boldepochs ; 
    end
    allboldtrials(:,1:16,:) = squeeze(allboldepochs(1,:,:,:)) ; allboldtrials(:,17:32,:) = squeeze(allboldepochs(2,:,:,:)) ; 
    colors = {'b','g','r'} ; 
    figure,
    for i=1:3
        shadedErrorBar([],mean(allboldtrials(i,:,:),2),std(allboldtrials(i,:,:),0,2)/sqrt(32),{colors{i}}) ; hold on ; 
    end
    title(subs{sub}) ; 
    allsubtrials(sub,:,:,:) = allboldtrials ; 
end

for i=1:3
   shadedErrorBar([],squeeze(mean(mean(allsubtrials(:,i,:,:),1),3)),squeeze(std(mean(allsubtrials(:,i,:,:),3),0,1))/sqrt(32),{colors{i}}) ; hold on ; xlim([1,21])
end
bsubtrials = allsubtrials - repmat(mean(allsubtrials(:,:,:,4:5),4),[1,1,1,size(allsubtrials,4)]) ; 
%mtspecs = squeeze(mean(allsubspecs(:,:,:,:,1:350),5)) ; 
bsubspecs = log(allsubspecs) - repmat(mean(log(allsubspecs(:,:,:,:,1:30)),5),[1,1,1,1,size(allsubspecs,5)]) ; 
mtbold = squeeze(mean(allsubtrials(:,:,:,12:14),4)) ;
%{
ztrials = zscore(squeeze(mean(allsubtrials(:,:,:,1:5),4)),[],3) ; 
% bad EEG trials
for i=1:8 ; 
    figure,
    for j=1:32 ; subplot(4,8,j) ; 
        imagesc(squeeze(bsubspecs(i,3,j,:,:)),[-2,2])
        title(j) ;
    end
end

% bad fmri trials
for i=1:8 
    figure
    for j=1:32
       subplot(4,8,j) ; 
       plot(squeeze(allsubtrials(i,3,j,:))) ; title(j) ; 
    end
end


ebads1 = {[19,21,27],[],[],[4,17,30],[22,24,30,32],[],[],[21,25,28]} ; 
ebads2 = {[10,15],[],[24],[13],[14,16,28],[],[],[3,5,6,17,26]} ; 
ebads3 = {[20],[],[15],[9,10,28],[12,17,18],[],[20],[7,29]} ; 
allebads = {ebads1,ebads2,ebads3} ; 
    
fbads1 = {[8,10,12,23,27],[10,26,30],[6,12,25],[4,9,17,29,31,26],[30],[19,28],[8,15,16,32],[19,27,32]} ; 
fbads2 = {[4,23],[8],[8,20,29],[31],[],[5,31],[1,29],[32]} ; 
fbads3 = {[4,11,30],[],[12,23,29,31],[],[26],[],[2,7,16,17,32],[20]} ; 
allfbads = {fbads1,fbads2,fbads2} ; 

for i=1:length(allebads)
    for j=1:length(allebads{i})
        ebads_ij = allebads{i}{j} ; 
        fbads_ij = allfbads{i}{j} ; 
        union_bads = union(ebads_ij,fbads_ij) ;
        goods = zeros(1,32) ; goods(union_bads) = 1 ; goodinds = find(goods==0) ; 
        unionbads{i,j} = goodinds ; 
    end
end


clear tcorrs
for i=1:8 ; 
    for j=1:3
        for k=1:63
            for el=1:size(allsubspecs,5)
                tcorrs(i,j,k,el) = corr2(squeeze(allsubspecs(i,j,unionbads{j,i},k,el)),squeeze(mtbold(i,j,unionbads{j,i}))) ; 
            end
        end
    end
end
%}
mtspecs = squeeze(mean(log(allsubspecs(:,:,:,:,50:280)),5)) ; 
mtbspecs = squeeze(mean((bsubspecs(:,:,:,:,50:280)),5)) ; 

clear mtcorrs
for i=1:8 ; 
    for j=1:3
        for k=1:63
            mtcorrs(i,j,k) = corr2(squeeze(mtspecs(i,j,:,k)),squeeze(mtbold(i,j,:))) ; 
        end
    end
end

colors = {'r','g','b'} ; figure,
for i=1:3
shadedErrorBar([],squeeze(mean(mtcorrs(:,i,:))),squeeze(std(mtcorrs(:,i,:),0,1))/sqrt(8),{colors{i}}) ; hold on ; xlim([1,63]) ; 
end
hline(0,'k') ; xlabel('frequency(hz)') ;ylabel('correlation (r)') ; 


eeg12(:,1:32,:) = squeeze(mtspecs(:,1,:,:)) ; eeg12(:,33:64,:) = squeeze(mtspecs(:,2,:,:)) ; 
eeg3 = squeeze(mtspecs(:,3,:,:)) ; 
fmri12(:,1:32) = squeeze(mtbold(:,1,:)) ; fmri12(:,33:64) = squeeze(mtbold(:,2,:)) ; 
fmri3 = squeeze(mtbold(:,3,:)) ; 
for i=1:8
    for j=1:size(eeg12,3)
        s12_corrs(i,j) = corr2(squeeze(eeg12(i,:,j)),fmri12(i,:)) ; 
        s3_corrs(i,j) = corr2(squeeze(eeg3(i,:,j)),fmri3(i,:)) ;         
    end
end
plot(mean(s12_corrs)) ; hold on ; plot(mean(s3_corrs),'r') ; 







