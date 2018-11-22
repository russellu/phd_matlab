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
    'tah'
    'thititip'
    'vincent'

} ; 
    
%%% load the raw FMRI images and extract the stimulus epochs
rawkey = 'common_*' ;
for sub=1:size(subjects,1) ; 
    cd(['c:/shared/allfmris/sub_',subjects{sub}]) ;
    ps = dir('pmat*') ; 
    size(ps) 
    mcorrs = load_nii('meancorrs.nii.gz') ; 
    cmask = mcorrs.img ; 
    post = zeros(size(cmask)) ; post(:,1:25,:) = 1 ; postmask = cmask.*post ; 
    occ = postmask > .25 ; 
    for i=1:size(ps,1) ; 
        nii = load_nii(ps(i).name) ; 
        pbrains{sub}(i,:,:,:) = nii.img ; 
    end
    occs{sub} = occ ; 
    
    
end


for i=1:size(pbrains,2)
    for j=1:size(pbrains{i},1) ; 
        c = 1-squeeze(pbrains{i}(j,:,:,:)) ; c(isnan(c)) = 0 ; 
        pmeans(i,j) = sum(sum(sum((c).*occs{i})))./sum(sum(sum(occs{i}))) ; 
        a = occs{i} ; 
        acount = 1 ; for ii=1:size(a,1) ; for jj=1:size(a,2) ; for k=1:size(a,3) ; if a(ii,jj,k)==1 ; spmats{i}(j,acount)  = c(ii,jj,k) ; acount = acount + 1 ; end ; end ; end ; end
    end
end
for i=1:size(spmats,2) ; meds(i,:) = median(spmats{i},2) ; end
bar(median(meds,1)) ; suptitle('median p-values for all subjects') ; 
set(gca,'XTick',1:6) ; set(gca,'XTickLabel',{'contrast','randomization','plaid','all','rnd lvls (1,3)'}) ; 


