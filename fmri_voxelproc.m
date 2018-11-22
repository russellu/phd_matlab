
clear all ; close all ; 
subs = {'sub_alex','sub_charest','sub_esteban','sub_fabio','sub_gab','sub_gabriella','sub_genevieve','sub_gina','sub_jeremie','sub_julie','sub_katrine','sub_lisa'...
        ,'sub_marc','sub_marie','sub_mathieu','sub_maxime','sub_mingham','sub_patricia','sub_po','sub_russell','sub_sunachakan','sub_vincent'} ;
for s=1:length(subs) ; 
    disp(subs{s}) ;
cd(['c:/shared/allfmris/',subs{s},'/ica']) ; ls 
denoised = load_untouch_nii('denoised.nii.gz') ; 
dimg = denoised.img ; 
clear epoched
icount = 1 ; 
for i=1:245:2205
    epoched(icount,:,:,:,:) = dimg(:,:,:,i:i+244) ;
    icount = icount + 1 ; 
end
hrf = spm_hrf(2) ; 
    cd ../trigs ; stimes = load('stimTimes1') ; stimes = stimes.stimTimes ;
    stimes = cell2mat(stimes) ; stimes = stimes(1:2:60) ; 
    trs = round(stimes/2) ; ts = zeros(1,245) ; ts(trs) = 1 ; 
    conved = conv(ts,hrf) ; conved = conved(1:245) ; 
    clear rvols ; 
    for i=1:size(epoched,1)
            rvols(i,:,:,:,:) = voxcorr(squeeze(epoched(i,:,:,:,:)),conved) ; 
    end
    stimes = dir('stimTimes*') ;
    clear alltimes alltrigs
    for st=1:length(stimes) ; 
        curr = load(stimes(st).name) ;
        curr = cell2mat(curr.stimTimes) ; 
        alltimes(st,:) = curr(1:2:60) ;
        alltrigs(st,:) = curr(2:2:60) ; 
    end
    stimcounts = [1,1,1,1,1,1] ; 
    alltimes = round(alltimes/2) ; 
    clear stimepochs
    for i=1:size(alltimes,1) ; 
        for j=1:size(alltimes,2)
            stimij = alltrigs(i,j) ; 
            epochij = squeeze(epoched(i,:,:,:,(alltimes(i,j)-1):(alltimes(i,j)+8))) ; 
            stimepochs(stimij,stimcounts(stimij),:,:,:,:) = epochij ; 
            stimcounts(stimij) = stimcounts(stimij) + 1 ; 
        end
    end   
    mtask = squeeze(mean(stimepochs(:,:,:,:,:,5:6),6)) ; 
    mrest = squeeze(mean(stimepochs(:,:,:,:,:,1:2),6)) ;     
    mdiff = squeeze(mean(mtask,2)-mean(mrest,2)) ; 
    f1 = load_untouch_nii('../f1.nii.gz') ;
    for i=1:6 ; 
        f1.img = squeeze(mdiff(i,:,:,:)) ; save_untouch_nii(f1,['vox_mdiff_',num2str(i),'.nii.gz']) ; 
        f1.img = squeeze(mean(mtask(i,:,:,:),2)) ; save_untouch_nii(f1,['vox_mtask_',num2str(i),'.nii.gz']) ; 
        f1.img = squeeze(mean(mrest(i,:,:,:),2)) ; save_untouch_nii(f1,['vox_mrest_',num2str(i),'.nii.gz']) ; 
    end 
    f1.img = squeeze(mean(rvols,1)) ; save_untouch_nii(f1,'rvols.nii.gz') ; 
end