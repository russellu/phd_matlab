%clear all ; close all 
%subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 
subs = {'MONG_01_RB','MONG_02_DP','MONG_03_CG','MONG_04_AB','MONG_05_SG','MONG_06_TS'} ; 

compinds = {[16,25],[8,23],[18,40],[4,21],[53,50],[17,20],[31,33]} ; 
badcomps = {[1,2,3,6,16,20,29,35:44,46:49,51:64],[1,2,4,5,7,10,11,16,28:64],[1,2,3,4,7,8,14,15,16,17,19,20,22,24,29,30,33,34:39,41:64],...
    [1,2,3,5,6,10:13,17:19,21:23,25:64],[1:6,9,10,14,16,20:25,27,28,30:33,35:38,41:64],[1:7,9,12,13,15,19:22,24,25,27:31,33,34,36:64],...
    [1:3,5,7,8,9,15,17,21,23,24:30,32,33,35:64]} ;

goodcomps = {[7,8,10,12,14,15,18,19,22,23,28,30],[3,8,12,13,18,25,26],[6,9,11,12,13,18,23,28,31,40],[7,8,20,24],[13,17,18,26,29,39],[8,17,18,23,32,35],[6,11,14,18,20,31,34]} ;
elecs = [59,45,31,46,60,9,20,18] ; 


for sub=1:length(subs)
    
    %%%% BOLD FMRI processing:
    %{
    cd(['c:/shared/badger_mri/',subs{sub},'/nii/melodic']) ; 
    mix = load('melodic_mix') ; 
    weights = load_untouch_nii('melodic_IC.nii.gz') ; 
    allweights{sub} = weights.img ; 
    meansub = load_untouch_nii('mean.nii.gz') ; 
    allmeans{sub} = meansub.img ; 
    mixinds = 1:735:735*5 ; 
    for i=1:length(mixinds) ; boldinds{i} = mixinds(i):mixinds(i)+734 ; end
    boldinds{6} = 735*5+1:735*5+1+449 ; 
    newboldinds{1} = boldinds{5} ; newboldinds{2} = boldinds{6} ; 
    boldinds = newboldinds ; 
    for i=1:length(boldinds) ; segmix{i} = mix(boldinds{i},:) ; end
    cd .. 
    %}
   % movie = load_untouch_nii('bp_*_retino_movie.nii.gz') ; 
   % rest = load_untouch_nii('bp_*_retino_rest.nii.gz') ; 

    %%%% EEG processing:
    cd(['c:\shared\mongoose\PROJET_MONGOOSE\raw\',subs{sub}]) ; 
    prefix = 'allfreq_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ;gammas2=dir([prefix,'preproc*gamma*set']) ;
    allstims1=dir([prefix,'preproc*allstim*set']) ;allstims2=dir([prefix,'preproc*allstim*set']) ;
    rests=dir('*REST*Pulse*vhdr') ; 
    movies=dir('*FIX*Pulse*vhdr') ;   
    %setNames = {allstims1(1).name,allstims1(2).name,gammas1(1).name,gammas1(2).name,movies(1).name,rests(1).name} ;   
    setNames = {movies(1).name,rests(1).name} ; 
    
    weights = pop_loadset('broad_merged.set') ; 
    hrf = spm_hrf(0.693) ; 
    
    for setN = 1:2;    
        EEG = pop_loadbv('.',setNames{setN}) ; 
        EEG = ica_applyweights(EEG,weights) ; 
        
       % bads = ones(1,64) ; bads(goodcomps{sub}) = 0 ; 
       % bads = find(bads==1) ; 
       % EEG = pop_subcomp(EEG,bads) ; 
        
        lats = {EEG.urevent.latency} ; 
        types = {EEG.urevent.type} ; 
        r128s = find(strcmp('R128',types)) ; 
        rlats = cell2mat(lats(r128s)) ; 
        
        eegacts = EEG.icaact(:,rlats(1):rlats(length(rlats))) ; 
        eegsecs = size(eegacts,2)./EEG.srate ; 
        eegTRs = round(eegsecs./0.693)+1 ; 
  
        [s,f] = spectopo(eegacts,0,EEG.srate,'plot','off') ; 
        
        eeghrf = spm_hrf(1/EEG.srate) ;
        
        filtcount = 1 ; clear filts
        for filt=1:2:100 ;
            filts(filtcount,:,:) = eegfiltfft(eegacts,EEG.srate,filt-2.5,filt+2.5) ;
            filtcount = filtcount + 1 ; 
        end

        clear convfilts allconvfilts 
        for i=1:64 ; disp(i) ; 
            for j=1:size(filts,1)
                convfilts = smooth(squeeze(filts(j,i,:)).^2,200) ; 
                convfilts = circshift(imresize(convfilts,[eegTRs,1]),0) ; 
                allconvfilts(i,j,:) = smooth(convfilts) ; 
            end
        end
        
        convpows{setN} = allconvfilts ; 
    
            meanfs(sub,setN,:,:) = s ; 

    end
    save('f','f') ; 
    %{
    
    % remove time points that are not in the EEG data:
    for i=1:length(segmix)
        segmix{i} = eegfiltfft(segmix{i}(1:size(convpows{i},3),:)',1/.693,.02,1)' ;
    end
    
    for i=1:length(segmix)
        for j=1:size(segmix{i},2)
            for k=1:size(convpows{i},1) 
                for el=1:size(convpows{i},2)
                    corrs(i,j,k,el) = corr2(squeeze(segmix{i}(40:end-40,j)),squeeze(convpows{i}(k,el,40:end-40))) ; 
                end
            end
        end
    end      
    
    figure,plot(squeeze(mat2gray(segmix{2}(40:end-40,25)))) ; hold on ; plot(squeeze(mat2gray(convpows{2}(10,9,40:end-40))),'r') ;
    %}
    restpow = convpows{2} ;  save('post_restpow','restpow') ; 
    moviepow = convpows{1} ; save('post_moviepow','moviepow') ; 
   
end
errorbar(squeeze(mean(mean(meanfs(:,1,elecs,1:200),1),3)),squeeze(std(mean(meanfs(:,1,elecs,1:200),3),0,1))./sqrt(7),'r') ; hold on ; 
errorbar(squeeze(mean(mean(meanfs(:,2,elecs,1:200),1),3)),squeeze(std(mean(meanfs(:,2,elecs,1:200),3),0,1))./sqrt(7),'b') ; 
set(gca,'XTick',1:20:200,'XTickLabel',round(f(1:20:200))) ; xlim([1,200]) ; xlabel('frequency(hz)') ; ylabel('log power') ; legend({'movie','rest'}) ; 

