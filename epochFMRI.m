%%%% MATLAB script to process FMRI data from visual experiments
% gets the stimulus times from log files, extracts stimulus type and time
% at which it occured, and saves an ideal convolved with an HRF to the
% subject specific directory
clear all ; close all

stimnames{1} = 'unperturbed'  ;
stimnames{2} = 'contrast_5%'  ;
stimnames{3} = 'contrast_33%'  ;
stimnames{4} = 'plaid'  ;
stimnames{5} = 'rnd_10%'  ;
stimnames{6} = 'rnd_60%'  ;

subjects = {
    'guillaume'
    'gina'
%    'maxime'
%    'marie'
%{
    'alex'
    'charest'
    'esteban'
    'fabio'
    'gab'
    'gabriella'
    'genevieve'
    'jeremie'
    'julie'
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
    'thititip'
    'vincent'
%}
} ; 

cd c:/shared/allfmris ; ls 
subs=dir('sub_*') ; 

hrf = [0,1,6,9,7,4,-1,0] ; 
ntrs = 245 ; 

for sub=1:size(subjects,1) ; % for all subjects
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
rawkey = 'band_common*' ;
for sub=1:size(subjects,1) ; 
    cd(['c:/shared/allfmris/sub_',subjects{sub}]) ;
    fmris = dir(rawkey) ; 
    fimgs = zeros(9,64,64,39,245) ; 
    corrnii = load_nii('postcorrs.nii.gz') ; 
    corrimg = corrnii.img ; 
    corrthresh = corrimg > .35 ; 
    negroi = load_nii('posroi.nii.gz') ;
    negroi = negroi.img ; 
    corrthresh = negroi ;
    
    clear fimgs
    for fmri=1:size(fmris,1)
        disp(fmris(fmri).name) ; 
        fnii = load_nii(fmris(fmri).name) ;   
        fimgs(fmri,:,:,:,:) = fnii.img ;  
    end
    clear subvoxs
    for fimg=1:size(fimgs,1) % for all 9 scans
        voxcount = 1 ; 
        for i=1:size(corrthresh,1)%x
            for j=1:size(corrthresh,2)%y
                for k=1:size(corrthresh,3)%z
                    if corrthresh(i,j,k)==1 % if this voxel was correlated with stimulus
                        subvoxs(fimg,voxcount,:) = fimgs(fimg,i,j,k,:) ; 
                        voxcount = voxcount + 1 ;
                    end
                end
            end
        end
    end
    meanvox = squeeze(mean(subvoxs,2)) ;
    sumvox = squeeze(sum(subvoxs,2)) ; 
    clear stimvals
    for scan=1:size(subvoxs,1)
         stimcounter = ones(1,6) ; 
        trigtimes = squeeze(alltrigs(sub,scan,:)) ; 
        trtimes = round(trigtimes/2) ; 
        for tr=1:size(trtimes,1)
            stimvals(scan,allparams(sub,scan,tr),stimcounter(allparams(sub,scan,tr)),:) = meanvox(scan,trtimes(tr)-1:trtimes(tr)+8) ; 
            stimcounter(allparams(sub,scan,tr)) = stimcounter(allparams(sub,scan,tr)) + 1 ; 
        end
    end
    % arrange all trials into a vector
    for stype=1:6 ; tcount = 1 ; for i=1:9 ; for j=1:5 ; alltrials(stype,tcount,:) = squeeze(stimvals(i,stype,j,:)) ; tcount = tcount + 1 ; end ; end  ; end 
    for i=1:size(alltrials,1)
        for j=1:size(alltrials,2)
            for k=1:size(alltrials,3)
                basenorms(i,j,k) = (alltrials(i,j,k) - alltrials(i,j,1)) ./ abs(alltrials(i,j,1)) ; 
            end        
        end
    end
    basenorms(isinf(basenorms)|isnan(basenorms)) = 0 ; 
    for i=1:6 ;
        [~,bads(i,:)] = iter_zthresh(squeeze(std(basenorms(i,:,:),0,3)),2) ;
    end
    clear goodtrials meantrials
    for i=1:size(basenorms,1)
        jcount = 1 ;
        for j=1:size(basenorms,2)
            if bads(i,j) == 0 
                goodtrials(jcount,:) = squeeze(basenorms(i,j,:)) ; 
                jcount = jcount + 1 ; 
            end     
        end
        allgoodtrials{sub,i} = goodtrials ; 
        meantrials(i,:) = mean(goodtrials,1) ; 
    end
    save('negrois','meantrials') ; 
    
    save('meantrials','meantrials') ; 
    
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

    
end

