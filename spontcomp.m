clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','sukhman','valerie'} ;
for subby=1%:length(subs) ; 
    cd(['c:/shared/badger_eeg/',subs{subby}]) ; %ls 
    % get the files in the proper order
    prefix = 'setrest_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ;gammas2=dir([prefix,'preproc*gamma*set']) ;
    allstims1=dir([prefix,'preproc*allstim*set']) ;allstims2=dir([prefix,'preproc*allstim*set']) ;
    movies=dir([prefix,'preproc*movie*set']) ; rests=dir([prefix,'preproc*rest*set']) ;
    %{
    if size(movies,1) ~= 0  
        setNames = {allstims1(1).name,allstims2(2).name,gammas1(1).name,gammas2(2).name,movies(1).name,rests(1).name} ; 
    else
        setNames = {allstims1(1).name,allstims2(2).name,gammas1(1).name,gammas2(2).name,rests(1).name} ; 
    end
    %}
        setNames = {rests(1).name} ; 

    for j=1:length(setNames) 
        eegs{j} = pop_loadset(setNames{j}) ; 
    end   
    % compute the frequency power
    clear corrs conversps mixscans
    for scan = 1 ; 
        clear ersp
        for c=1:64 
            dummies = round(0.693*3) ; %3*TR dummies
            ntrs = ceil(eegs{scan}.xmax/0.693) ; 
            [ersp(c,:,:),itc,powbase,times,freqs,~,~] = newtimef(eegs{scan}.icaact(c,:),eegs{scan}.pnts,[eegs{scan}.xmin,eegs{scan}.xmax],eegs{scan}.srate,0,...
                'plotersp','off','plotitc','off','baseline',NaN,'timesout',ntrs,'nfreqs',120,'freqs',[1,120],'winsize',200) ; 
        end
        hrf = spm_hrf(0.693) ; 
        conversp = zeros(size(ersp,1),size(ersp,2),size(ersp,3)) ; 
        for c=1:size(ersp,1)
            for f=1:size(ersp,2)
                 conv_cf = conv(squeeze(ersp(c,f,:)),hrf,'full') ; 
                 conv_cf = conv_cf(1:size(ersp,3)) ; 
                 conversp(c,f,:) = conv_cf ; 
            end
        end
        % get the FMRI
        %cd(['c:/shared/badger_mri/',subs{subby},'/nii/melodic']) ;
        
         cd(['c:/shared/badger_mri/',subs{subby},'/nii/melodic']) ;
         ls
        indends = 1:734:735*5;
        inds = {1:indends(2),indends(2):indends(3),indends(3):indends(4),indends(4):indends(5),indends(5):indends(6),indends(6):indends(6)+449} ;
        mix = load('melodic_mix') ;
        mixes{scan} = mix(inds{scan},:) ;  
        mixscan = mixes{scan}(1:size(ersp,3),:) ; 
        mixscan = eegfiltfft(mixscan',1/0.693,0.02,1)' ; 
        for f=1:size(ersp,2)
            corrs(scan,f,:,:) = corr(squeeze(conversp(:,f,15:end-15))',mixscan(15:end-15,:)) ; 
        end  
        conversps{scan} = conversp ; 
        mixscans{scan} = mixscan ; 
    end   
    %for i=1:size(corrs,4) ; figure,imagesc(squeeze(mean(corrs(:,:,:,i))),[-0.5,.5]) ; end
    %comps = [14,23,24] ; 
    %{
    for i=1:size(corrs,4) ; figure,
        for j=1:6 ; subplot(2,3,j) ; 
            imagesc(imfilter(squeeze(mean(corrs(j,:,:,i),1)),fspecial('gaussian',[3,1],1.5)),[-0.75,0.75]) ;
        end
    end
    %figure,  %  plot(squeeze(mean(mean(corrs(1:2,:,comps,:),1),3)),'r') ; hold on ;  plot(squeeze(mean(mean(corrs(3:4,:,comps,:),1),3)),'k') ;
    %for i=1:3 ; subplot(2,2,i) ; imagesc(squeeze(mean(corrs(1:2,:,comps(i),:))),[-.25,.25]) ; end
    eegcomp = 18;ff=36;fmricomp=48 ; tcount=2; 
    
    
    figure,plot(mat2gray(squeeze(conversps{tcount}(eegcomp,ff,15:end-15)))) ; hold on ; plot(mat2gray(mixscans{tcount}(15:end-15,fmricomp)),'r');
    figure,plot(mat2gray(squeeze(conversps{tcount}(eegcomp,ff,15:end-15))),mat2gray(mixscans{tcount}(15:end-15,fmricomp)),'o'); 
    title(num2str(corr(mat2gray(squeeze(conversps{tcount}(eegcomp,ff,15:end-15))),mat2gray(mixscans{tcount}(15:end-15,fmricomp)))))
    %}
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(corrs(:,:,i,:),1)),[-.5,.5]) ; end

    %[c,p] = corr(mat2gray(squeeze(conversp(eegcomp,ff,15:end-15))),mat2gray(mix1(15:end-15,fmricomp))) ; 
    %setNames{6}
    %EEG = pop_loadset(preprocs(i).name) ; 
    %if i==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged,1) ; end 
end

