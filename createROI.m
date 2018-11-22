clear all ; close all ;
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
for s=8
    cd(['c:/shared/badger_eeg2/',subs{s}]) ;
    bcgs = dir('bcgica*gamma*set') ;

    for b=1:length(bcgs)
        %if b==1 ; name = dir('remove_bcg*gamma*1*set') ; else name = dir('remove_bcg*gamma*2*set') ; end
        EEG = pop_loadset(bcgs(b).name) ;     
        EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;

        triglabs = {'S  1','S  2','S  3'} ; 
        lats = cell2mat({EEG.urevent.latency}) ; 
        trigs = {EEG.urevent.type} ; r128s = find(strcmp('R128',trigs)) ;
        startlat = lats(r128s(1)) ; 
        for t=1:length(triglabs)
            alltrigs = find(strcmp(triglabs{t},trigs)) ; 
            allLats(t,:) = lats(alltrigs) ; 
        end
        allLats = allLats - startlat ; 
        roundres = round(allLats/(0.693*EEG.srate)) ; roundres = roundres(:) ;
        zts = zeros(1,735) ; 
        for r=1:length(roundres) ; zts(roundres(r):roundres(r)+round(5/0.693)) = 1 ; end
        hrf = spm_hrf(0.693) ; 
        conved = conv(zts,hrf,'full') ; conved = conved(1:735) ; 
        allconved(b,:) = conved ; 
    end
    
    cd(['c:/shared/newbadger_mri/',subs{s}]) ; 
    bps = dir('bp_*gamma*') ; 
    for bp=1:length(bps)
        nii = load_untouch_nii(bps(bp).name) ; 
        corrs(:,:,:,bp) = voxcorr(nii.img(:,:,:,40:end-40),allconved(bp,40:end-40)) ; 
    end
    f1=  load_untouch_nii('f1.nii.gz') ; f1.img = mean(corrs,4) ; 
    save_untouch_nii(f1,'corrs.nii.gz') ; 
end