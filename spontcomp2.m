clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','valerie'} ;
for subby=1%:length(subs) ; 
    cd(['c:/shared/badger_eeg/',subs{subby}]) ; %ls 
    % get the files in the proper order
    prefix = '' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ;gammas2=dir([prefix,'preproc*gamma*set']) ;
    allstims1=dir([prefix,'preproc*allstim*set']) ;allstims2=dir([prefix,'preproc*allstim*set']) ;
    movies=dir([prefix,'preproc*movie*set']) ; rests=dir([prefix,'preproc*rest*set']) ;

    setNames = {gammas1(1).name,gammas2(1).name,allstims1(1).name,allstims2(1).name} ; 
    fNames = {'reg_mc_retino_gamma_01.nii.gz','reg_mc_retino_gamma_02.nii.gz','reg_mc_retino_allstims_01.nii.gz','reg_mc_retino_allstims_02.nii.gz'} ; 
    for j=1:length(setNames) 
        eegs{j} = pop_loadset(setNames{j}) ; 
    end   
    % compute the frequency power
    clear corrs conversps mixscans
    for scan = 1:4 ; 
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
        

        cd(['c:/shared/badger_mri/',subs{subby},'/nii/']) ;
        rest = load_untouch_nii(fNames{scan}) ; 
        rest.img = rest.img(:,:,:,1:size(conversp,3));
        size(rest)
        for f=1:64 
           corrs(scan,:,:,:,f) = voxcorr(rest.img(:,:,:,50:end-15),squeeze(mean(conversp(f,8:25,50:end-15),2))) ;          
        end
     
    end   
    anat = load_untouch_nii('fwarp.nii.gz') ; 
    mcorrs = squeeze(mean(corrs(1:4,:,:,:,:),1)) ;
    %figure,for i=1:64 ; subplot(5,13,i) ; imagesc(imrotate(squeeze(mean(mcorrs(:,:,15:25,i),3)),90),[-.25,.25]) ; title(i) ; end
    zinds = 23:28 ; figure,
    for x=1:size(corrs,4) ; subplot(5,13,x),
        plotoverlayIntensity2D(squeeze(mean(anat.img(zinds,:,1:end-2),1)),mat2gray(abs(squeeze(mean(mcorrs(zinds,:,1:end-2,x),1)))),squeeze(mean(mcorrs(zinds,:,1:end-2,x),1)),90) ;  title(x) ; 
    end
    suptitle(subs{subby}); 
end

