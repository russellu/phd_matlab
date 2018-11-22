clear all ; close all

subs = {'alex','dina','genevieve','jeremie','karl','russell','valerie','tegan'} ; 
prefix='allfreq_' ; 

for sb=1:length(subs) ; 
    cd(['c:/shared/badger_eeg/',subs{sb}]) ; ls 
    gammas=dir([prefix,'preproc*gamma*set']) ;
    for g=1:2 
        clear stims
        cd(['c:/shared/badger_eeg/',subs{sb}]) ; ls 
        EEG = pop_loadset(gammas(g).name) ; 
        
        trigsr = {EEG.urevent.type} ; latsr = cell2mat({EEG.urevent.latency}) ; 
        r128s = find(strcmp(trigsr,'R128')) ;
        rlats = latsr(r128s) ; 
        EEG = pop_select(EEG,'nopoint',[1,rlats(1) ; rlats(length(rlats)),size(EEG.data,2)]) ; 
        
        events = {EEG.urevent.type} ; 
        lats = {EEG.urevent.latency} ;
        stims{1} = find(strcmp('S  1',events)) ; stims{2} = find(strcmp('S  2',events)) ; stims{3} = find(strcmp('S  3',events)) ; 
        allinds = [cell2mat(stims(1)),cell2mat(stims(2)),cell2mat(stims(3))] ; 
        [~,si] = sort(allinds) ; 
        stimlats = cell2mat(lats(allinds(si))) ; 

        fmris = find(strcmp('R128',events)) ; 
        fmrilats = cell2mat(lats(fmris)) ; 
        
        rawstimlats = stimlats - fmrilats(1) ; 

        % find the fmri TR indices to which each stimulus is closest:
        design = zeros(1,length(fmris)) ; 
        longdesign = zeros(1,size(EEG.data,2)+round(EEG.srate/0.693)) ; 
        for i=1:length(stimlats)
           longdesign(rawstimlats(i):rawstimlats(i)+5*EEG.srate) = 1 ; 
           lati = find(abs(fmrilats-stimlats(i))==min(abs(fmrilats-stimlats(i)))) ;  
           design(lati:lati+round(5/0.693)) = 1 ;
           vlines(i,:) = lati:lati+round(5/0.693) ; 
        end
        shortdesign = imresize(longdesign,[1,length(fmris)]) ; 
        if length(shortdesign) < 735 ; shortdesign(length(shortdesign)+1:735) = mean(shortdesign) ; end % if you somehow closed the EEG before the FMRI was done. lol
        cd(['c:/shared/badger_mri/',subs{sb},'/nii',]) ; ls 
        
       % nii = load_untouch_nii('bp_reg_topup_mc_retino_gamma_01.nii.gz') ; 
        hrf = spm_hrf(0.693) ; 
        convshort = conv(shortdesign,hrf,'full') ; convshort = convshort(1:length(shortdesign)) ;
        dlmwrite(['gamma_design_',num2str(g),'.txt'],shortdesign') ; 
        dlmwrite(['gamma_hrf_',num2str(g),'.txt'],convshort') ; 
       % vline(vlines(:)) ; hold on ; plot(mat2gray(a(:,64)),'LineWidth',2) ;  plot(mat2gray(convshort(:)),'r','LineWidth',2) ; xlabel('tr') ; xlim([35,690])
       % title(['r=',num2str(corr(a(35:690,64),convshort(35:690)'))]) ;    
    end
    
end











