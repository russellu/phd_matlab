%%%% MATLAB script to process FMRI data from visual experiments
% gets the stimulus times from log files, extracts stimulus type and time
% at which it occured, and saves an ideal convolved with an HRF to the
% subject specific directory
clear all ; close all

stimnames{1} = 'unperturbed'  ;
stimnames{2} = 'contrast=5%'  ;
stimnames{3} = 'contrast=33%'  ;
stimnames{4} = 'plaid'  ;
stimnames{5} = 'rnd=10%'  ;
stimnames{6} = 'rnd=60%'  ;

subjects = {
    'alex'
    'charest'
    'esteban'
    'fabio'
    'gab'
    'gabriella'
    'genevieve'
    'gina'
    'guillaume'
    'jeremie'
    'julie'
    'katrine'
    'lisa'
    'marc'
    'marie'
    'mathieu'
    'mingham'
    'maxime'
    'patricia'
    'po'
    'russell'
    'sunachakan'
    'tah'
    'thititip'
    'vincent'
} ; 

cd c:/shared/allfmris ; ls 
subs=dir('sub_*') ; 

hrf = spm_hrf(2) ; %[0,1,6,9,7,4,-1,0] ; 
ntrs = 245 ; 

for sub=23%:size(subjects,1) ; % for all subjects
    cd(['c:/shared/allfmris/sub_',subjects{sub}]) ;
    cd trigs ; 
    alltimes = dir('stimTimes*.mat') ; 
    for stimFile=1:size(alltimes,1) ; % for all scans (1-9 usually)
        ctime = load(alltimes(stimFile).name) ;
        scanTimes = ctime.stimTimes ; 
        for presentation=1:size(scanTimes,2) % for all presentations within the scan
            alltrigs(sub,stimFile,presentation) = scanTimes{presentation}(1) ; 
            allparams(sub,stimFile,presentation) = scanTimes{presentation}(2) ; 
        end
    end
    % create ideal from first time series and convolve with canonical HRF
    stimes = round(squeeze(alltrigs(sub,1,:))) ; % trigger presentation times
    trtimes = stimes/2 ; % converted to TR
    ideal = zeros(1,ntrs) ; 
    for i=1:size(trtimes,1) % fill the ideal with 1s where there was a stimulus
        ideal(trtimes(i)) = 1 ; 
    end
    idealhrf = conv(ideal,hrf) ; % convolve ideal with HRF
    idealhrf = idealhrf(1:245) ; % shave off last few indices
    cd .. ; dlmwrite('conved',idealhrf') ; % save in subject's directory
end
        
%%% load the raw FMRI images and extract the stimulus epochs
rawkey = 'common_*' ;
for sub=23%size(subjects,1) ; 
    cd(['c:/shared/allfmris/sub_',subjects{sub}]) ;
    fmris = dir(rawkey) ; 
    fimgs = zeros(9,64,64,39,245) ; 
    %corrnii = load_nii('postcorrs.nii.gz') ; 
    %corrimg = corrnii.img ; 
    %corrthresh = corrimg > .3 ;
    %posroi = load_nii('posroi.nii.gz') ;
    %posroi = posroi.img ; 
    %negroi = load_nii('negroi.nii.gz') ;
    %negroi = negroi.img ; 
    %posroi = corrimg > .35 ; 
    %negroi = corrimg < -.15 ; 
    
    clear fimgs
    for fmri=1:size(fmris,1)
        disp(fmris(fmri).name) ; 
        fnii = load_nii(fmris(fmri).name) ;   
        fimgs(fmri,:,:,:,:) = fnii.img ;  
    end
    
    
    
    mcorrs = load_nii('meancorrs.nii.gz')  ;
    cmask = mcorrs.img ;  
    s = zeros(size(cmask)) ; s(:,1:30,:) = 1 ; cmask = cmask.*s ;

    
    
    %{
    brainepochs = zeros(6,size(fimgs,2),size(fimgs,3),size(fimgs,4),45,10) ; 
    subparams = squeeze(allparams(sub,:,:)) ; 
    stimcounts = ones(1,6) ;
    for scan=1:size(subparams,1) 
        for stim=1:size(trtimes,1) ; 
            brainepochs(subparams(scan,stim),:,:,:,stimcounts(subparams(scan,stim)),:) = ...
                                            fimgs(scan,:,:,:,(trtimes(stim)-2):trtimes(stim)+7) ; 
            stimcounts(subparams(scan,stim)) = stimcounts(subparams(scan,stim)) + 1 ; 
        end
    end
    clear percbrain
    %for i=1:6 ; for j=1:45 ; ; percbrain{i}(:,:,:,j) = (squeeze(mean(brainepochs(i,:,:,:,j,3:7),6))-squeeze(mean(brainepochs(i,:,:,:,j,1:2),6)))./squeeze(mean(brainepochs(i,:,:,:,j,1:2),6)) ; end ; end
    for i=1:6 ; for j=1:45 ;  percbrain(i,:,:,:,j) = (squeeze(mean(brainepochs(i,:,:,:,j,3:7),6))-squeeze(mean(brainepochs(i,:,:,:,j,1:2),6)))./squeeze(mean(brainepochs(i,:,:,:,j,1:2),6)) ; end ; end
    mpercbrain = squeeze(mean(percbrain,5)) ; mpercbrain(isnan(mpercbrain)) = 0 ; mpercbrain(isinf(mpercbrain)) = 0 ; mpercbrain(mpercbrain > .2) = 0 ; 
    for i=1:6 ; mcorrs.img = squeeze(mpercbrain(i,:,:,:)) ; save_nii(mcorrs,['perc_',num2str(i),'.nii.gz']) ; end
    %}
    
 %   clear mepochs ; for i=1:6 ; mepochs(i,:,:,:,:) = mean(squeeze(brainepochs(i,:,:,:,:,:)-repmat(mean(brainepochs(i,:,:,:,:,:),6),[1,1,1,1,1,10])),4) ; end
    
%%%% remove the bad trials %%%%
%{
    mcorrs = load_nii('meancorrs.nii.gz')  ;
    cmask = mcorrs.img ; 
    post = zeros(size(cmask)) ; post(:,1:25,:) = 1 ; postmask = cmask.*post ; 
    occ = postmask > .3 ; 
    clear allvox
    for i=1:size(percbrain,2) ;
        for j=1:size(percbrain{i},4) ;
            curbrain = percbrain{i}(:,:,:,j) ; 
            allvox(i,j,:) = curbrain(find(occ==1)) ; 
        end
    end
    xs = (squeeze(sum(diff(allvox(:,:,:),1,3).^2,3))) ;
    for i=1:6 ; [~,bads(i,:)] = iter_zthresh(xs(i,:),3) ; end
    for i=1:6 ;
        s = squeeze(mean(percbrain{i}(:,:,:,find(bads(i,:)==0)),4)) ;
        holder = zeros(size(s,1),size(s,2),size(s,3),1) ; holder(:,:,:,1) = s ; holder = repmat(holder,[1,1,1,sum(bads(i,:))]) ; 
        percbrain{i}(:,:,:,find(bads(i,:)==1)) = holder ;
    end % interpolate bad trials with mean of good trials

%}
    %percbrain = squeeze((mean(brainepochs(:,:,:,:,:,4:8),6)-mean(brainepochs(:,:,:,:,:,1:2),6))./mean(brainepochs(:,:,:,:,:,1:2),6)) ; 
    %{
    stimtypes = {[1,3,2],[1,5,6],[1,4],[1:6],[1,6]} ; 
    clear pmat
   % for s=1:size(stimtypes,2)
        for i=1:size(brainepochs,2)
            for j=1:size(brainepochs,3)
                for k=1:size(brainepochs,4)
                    %{
                    pmat(1,i,j,k) = anova1([squeeze(percbrain{1}(i,j,k,:)),squeeze(percbrain{3}(i,j,k,:)),squeeze(percbrain{2}(i,j,k,:))] ,[],'off') ; 
                    pmat(2,i,j,k) = anova1([squeeze(percbrain{1}(i,j,k,:)),squeeze(percbrain{5}(i,j,k,:)),squeeze(percbrain{6}(i,j,k,:))] ,[],'off') ; 
                    pmat(3,i,j,k) = anova1([squeeze(percbrain{1}(i,j,k,:)),squeeze(percbrain{4}(i,j,k,:))] ,[],'off') ; 
                    pmat(4,i,j,k) = anova1([squeeze(percbrain{1}(i,j,k,:)),squeeze(percbrain{2}(i,j,k,:)),squeeze(percbrain{3}(i,j,k,:)),squeeze(percbrain{4}(i,j,k,:)),squeeze(percbrain{5}(i,j,k,:)),squeeze(percbrain{6}(i,j,k,:))] ,[],'off') ; 
                    pmat(5,i,j,k) = anova1([squeeze(percbrain{1}(i,j,k,:)),squeeze(percbrain{6}(i,j,k,:))] ,[],'off') ; 
                    %}
                    [~,f] = anova1([squeeze(percbrain(1,i,j,k,:)),squeeze(percbrain(3,i,j,k,:)),squeeze(percbrain(2,i,j,k,:))] ,[],'off') ; 
                    if ~isempty(f{2,5}) ; pmat(1,i,j,k) = f{2,5} ; else pmat(1,i,j,k) = 0 ; end
                    [~,f] = anova1([squeeze(percbrain(1,i,j,k,:)),squeeze(percbrain(5,i,j,k,:)),squeeze(percbrain(6,i,j,k,:))] ,[],'off') ; 
                    if ~isempty(f{2,5}) ; pmat(2,i,j,k) = f{2,5} ; else pmat(2,i,j,k) = 0 ; end
                    [~,f] = anova1([squeeze(percbrain(1,i,j,k,:)),squeeze(percbrain(4,i,j,k,:))] ,[],'off') ; 
                    if ~isempty(f{2,5}) ; pmat(3,i,j,k) = f{2,5} ; else pmat(3,i,j,k) = 0 ; end
                    [~,f] = anova1([squeeze(percbrain(1,i,j,k,:)),squeeze(percbrain(2,i,j,k,:)),squeeze(percbrain(3,i,j,k,:)),squeeze(percbrain(4,i,j,k,:)),squeeze(percbrain(5,i,j,k,:)),squeeze(percbrain(6,i,j,k,:))] ,[],'off') ;  
                    if ~isempty(f{2,5}) ; pmat(4,i,j,k) = f{2,5} ; else pmat(4,i,j,k) = 0; end
                    [~,f] = anova1([squeeze(percbrain(1,i,j,k,:)),squeeze(percbrain(6,i,j,k,:))] ,[],'off') ; 
                    if ~isempty(f{2,5}) ; pmat(5,i,j,k) = f{2,5} ; else pmat(5,i,j,k) = 0 ; end
                end
            end
        end
    %end
    
    
    
   
   a = {squeeze(percbrain{1}(20,20,20,:)),squeeze(percbrain{1}(20,20,20,:)),squeeze(percbrain{1}(20,20,20,:))} ; 
    
    for i=1:size(pmat,1)
        mcorrs.img = squeeze(pmat(i,:,:,:)) ; 
        save_nii(mcorrs,['pmat_',num2str(i),'.nii.gz']) ; 
    end
   
    %}
    
    disp(subjects{sub}) ; 
   % pm = squeeze(mean(percbrain,5)) ; 
   % for i=1:6 ; for j=1:45 ; percs(i,j) = squeeze(sum(sum(sum(squeeze(percbrain(i,:,:,:,j)).*cmask)))) ; end ; end
   % mcorrs.img = squeeze(1-pmat) ;
   % save_nii(mcorrs,'pbrain.nii.gz') ; 

    thresh = 0.35 ; 
    roi = cmask ; %load_nii('occroi.nii.gz') ; roi = roi.img ; 
    clear negsubvoxs possubvoxs
    for fimg=1:size(fimgs,1) % for all 9 scans
        posvoxcount = 1 ; negvoxcount = 1 ;
        for i=1:size(cmask,1)%x
            for j=1:size(cmask,2)%y
                for k=1:size(cmask,3)%z
                   % if negroi(i,j,k)==1 && roi(i,j,k) ==1 % if this voxel was part of the negative ROI
                   %     negsubvoxs(fimg,negvoxcount,:) = fimgs(fimg,i,j,k,:) ; 
                   %     negvoxcount = negvoxcount + 1 ;
                   % end
                    if cmask(i,j,k)>=thresh %&& roi(i,j,k) ==1 % if this voxel was part of the positive ROI
                        possubvoxs(fimg,posvoxcount,:) = fimgs(fimg,i,j,k,:) ; 
                        posvoxcount = posvoxcount + 1 ;
                    end
                end
            end
        end
    end
    
    subparams = squeeze(allparams(sub,:,:)) ; 
    posepochs = zeros(6,size(possubvoxs,2),45,11) ; 
    %negepochs = zeros(6,size(negsubvoxs,2),45,11) ; 
    
    stimcounts = ones(1,6) ;
    for scan=1:9
        for trial=1:30
            posepochs(subparams(scan,trial),:,stimcounts(subparams(scan,trial)),:) = possubvoxs(scan,:,(trtimes(trial)-2):trtimes(trial)+8) ; 
            stimcounts(subparams(scan,trial)) = stimcounts(subparams(scan,trial)) + 1 ; 
        end   
    end
 %   stimcounts = ones(1,6) ;
 %   for scan=1:9
 %       for trial=1:30
 %           negepochs(subparams(scan,trial),:,stimcounts(subparams(scan,trial)),:) = negsubvoxs(scan,:,(trtimes(trial)-1):trtimes(trial)+9) ; 
 %           stimcounts(subparams(scan,trial)) = stimcounts(subparams(scan,trial)) + 1 ; 
 %       end   
 %   end
    save('posepochs','posepochs') ; 
    
    
    
   % allpos{sub}(:,:,:,:) = posepochs ;
 %   allneg{sub}(:,:,:,:) = negepochs ; 

    %{
    
    %stimspecfmri.m requires that epoch_fmri is run previous.
    % gets the stimulus specific activations maps. 
    % needs the variable fimgs to be set
    % pooledvar = (var1*(n1-1) + var2*(n2-1)) / (n1+n2-2) 
    % SE = pooledvar * sqrt(1/n1 + 1/n2) ;
    clear allstims ; 
    tcounts = ones(6,1); 
    for trial=1:size(allparams,2) ;
        for s=1:size(trtimes,1)
            allstims(allparams(sub,trial,s),tcounts(allparams(sub,trial,s)),:,:,:,:) = fimgs(trial,:,:,:,trtimes(s)-2:trtimes(s)+4) ; 
            tcounts(allparams(sub,trial,s)) = tcounts(allparams(sub,trial,s)) + 1 ;
        end
    end
    tvals = (squeeze(mean(mean(allstims(:,:,:,:,:,5:7),2),6)-mean(mean(allstims(:,:,:,:,:,1:2),2),6)))./squeeze(std(mean(allstims(:,:,:,:,:,3:7),6),0,2)./sqrt(45)) ;
    save('tvals','tvals') ;
    baset = squeeze(mean(mean(allstims(:,:,:,:,:,1:2),2),6)) ; 
    taskt = squeeze(mean(mean(allstims(:,:,:,:,:,3:7),2),6)) ; 
    
    varbase = squeeze(var(mean(allstims(:,:,:,:,:,1:2),6),0,2)) ; 
    vartask = squeeze(var(mean(allstims(:,:,:,:,:,3:7),6),0,2)) ; 
    pooledvar = (varbase*44 + vartask*44) ./ 88 ; 
    pooledstd = sqrt(pooledvar) ; 
    SE = pooledstd.* sqrt(1/44 + 1/44) ; 
    t = (taskt - baset) ./ SE ; 
    
    x = load_nii('meancorrs.nii.gz') ; 
    
    for i=1:6 ; 
        x.img = squeeze(t(i,:,:,:)) ; 
        save_nii(x,['stim_',stimnames{i},'.nii.gz']) ;
    end
    
    v = squeeze(mean(varbase,1)) ; x.img = v ; save_nii(x,'varbase.nii.gz') ; 
    v = squeeze(mean(vartask,1)) ; x.img = v ; save_nii(x,'vartask.nii.gz') ; 
    v = squeeze(mean(varbase-vartask,1)) ; x.img = v ; save_nii(x,'varbase_min_vartask.nii.gz') ; 
    %}
end

%stypes = squeeze(mean(allpos{1},2)) ;
%bcorrs = (stypes(:,:,:) - repmat(mean(stypes(:,:,1:2),3),[1,1,11]))./repmat(mean(stypes(:,:,1:2),3),[1,1,11]) ; for i=1:6 ; subplot(2,3,i) ; imagesc(squeeze(bcorrs(i,:,:)),[-.02,.02]) ; end
%errorbar(squeeze(mean(bcorrs,2))',squeeze(std(bcorrs,0,2)./sqrt(45))','LineWidth',2) ; legend(stimnames) ;
%{
for i=1:size(allpos,2)
    posvoxs(i,:,:,:) = squeeze(mean(allpos{i},2)) ; 
    negvoxs(i,:,:,:) = squeeze(mean(allneg{i},2)) ; 
end
% normalize the negative responses
for i=1:size(allpos,2) ; 
    for j=1:6 ; 
        for k=1:45 ; 
            normpos(i,j,k,:) = (posvoxs(i,j,k,:)-squeeze(mean(posvoxs(i,j,k,1),4)))./squeeze(mean(posvoxs(i,j,k,1),4)) ; 
            normneg(i,j,k,:) = (negvoxs(i,j,k,:)-squeeze(mean(negvoxs(i,j,k,1),4)))./squeeze(mean(posvoxs(i,j,k,1),4)) ; 
        end ; 
    end ;
end

for i=1:size(subjects,1) ; 
    cd(['c:/shared/allfmris/sub_',subjects{i}]) ;
    npi = squeeze(normpos(i,:,:,:)) ;
    nni = squeeze(normneg(i,:,:,:)) ; 
    save('npi','npi') ; 
    save('nni','nni') ; 
end
%}
