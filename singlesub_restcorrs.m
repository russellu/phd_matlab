clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 

for sub=1%:length(subs)
    
    %%%% BOLD FMRI processing:
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
    cd .. 
    movie = load_untouch_nii('bp_*_retino_movie.nii.gz') ; 
    rest = load_untouch_nii('bp_*_retino_rest.nii.gz') ; 

    %%%% EEG processing:
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    prefix = 'allfreq_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ;gammas2=dir([prefix,'preproc*gamma*set']) ;
    allstims1=dir([prefix,'preproc*allstim*set']) ;allstims2=dir([prefix,'preproc*allstim*set']) ;
    movies=dir([prefix,'preproc*movie*set']) ; rests=dir([prefix,'preproc*rest*set']) ;
    setNames = {allstims1(1).name,allstims1(2).name,gammas1(1).name,gammas1(2).name,movies(1).name,rests(1).name} ;   

    hrf = spm_hrf(0.693) ; 
    
    for setN = 1:6;    
        EEG = pop_loadset(setNames{setN}) ; 
        lats = {EEG.urevent.latency} ; 
        types = {EEG.urevent.type} ; 
        r128s = find(strcmp('R128',types)) ; 
        rlats = cell2mat(lats(r128s)) ; 
        
        eegacts = EEG.icaact(:,rlats(1):rlats(length(rlats))) ; 
        eegsecs = size(eegacts,2)./EEG.srate ; 
        eegTRs = round(eegsecs./0.693)+1 ; 
    
        %{
        clear epow
        for c=1:64%length(comps{sub})  
            [epow(c,:,:),itc,powbase,times,freqs,~,~] = newtimef(eegacts(c,:),size(eegacts,2),[0,size(eegacts,2)./EEG.srate],EEG.srate,0,...
                'plotersp','off','plotitc','off','baseline',NaN,'timesout',eegTRs,'nfreqs',50,'winsize',250,'freqs',[1,100]) ; 
        end
        allpows{setN} = epow ; 
        %}
        
        eeghrf = spm_hrf(1/EEG.srate) ;
        
        low = eegfiltfft(eegacts,EEG.srate,1,7) ;
        alphabeta = eegfiltfft(eegacts,EEG.srate,8,25) ; 
        gamma = eegfiltfft(eegacts,EEG.srate,30,80) ;      
        
        clear allconvalphabeta allconvgamma allconvlow
        for i=1:64
            %convalphabeta = conv(abs(alphabeta(i,:)),eeghrf,'full') ; convalphabeta = convalphabeta(1:length(alphabeta)) ; 
            convalphabeta = smooth(alphabeta(i,:).^2,200) ; 
            convalphabeta = circshift(imresize(convalphabeta,[eegTRs,1]),6) ; 
            allconvalphabeta(i,:) = smooth(convalphabeta) ; 
            convgamma = smooth(gamma(i,:).^2,200) ; 
            convgamma = circshift(imresize(convgamma,[eegTRs,1]),6) ; 
            allconvgamma(i,:) = smooth(convgamma) ; 
            convlow = smooth(low(i,:).^2,200) ; 
            convlow = circshift(imresize(convlow,[eegTRs,1]),6) ; 
            allconvlow(i,:) = smooth(convlow) ; 
        end
        
        convpows{setN}(1,:,:) = allconvlow ; 
        convpows{setN}(2,:,:) = allconvalphabeta ; 
        convpows{setN}(3,:,:) = allconvgamma ; 
    end
    
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
    
    figure,plot(squeeze(mat2gray(segmix{6}(40:end-40,67)))) ; hold on ; plot(squeeze(mat2gray(convpows{6}(2,32,40:end-40))),'r') ;
    
    %plot(squeeze((segmix{4}(40:end-40,64))),smooth(squeeze((convpows{3}(3,30,40:end-40)))),'.') ; title(corr2(squeeze((segmix{3}(40:end-40,64))),smooth(squeeze((convpows{4}(3,30,40:end-40))))))
    %{
    hrf = spm_hrf(0.693) ; 
    for i=1:length(allpows)
       for j=1:size(allpows{i},1)
           for k=1:size(allpows{i},2)
                convijk = conv(squeeze(allpows{i}(j,k,:)),hrf,'full') ; 
                convpows{i}(j,k,:) = convijk(1:size(allpows{i},3)) ; 
           end
       end
    end
    for i=1:6 
        segmixi = segmix{i} ; 
        epowi = allpows{i} ; 
        for j=1:size(segmixi,2)
            for k=1:size(epowi,1)
                for el=1:50
                    corrs(i,j,k,el) = corr2(segmix{i}(40:end-40,j),squeeze(convpows{i}(k,el,40:end-40))) ;
                end
            end
        end
    end
      %}
    
end

%figure,plot(squeeze(mean(convpows{3}(26,25:30,:),2))+squeeze(mean(convpows{3}(30,25:30,:),2)))

