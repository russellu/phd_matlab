clear all ; close all 
subs = {'sub_alex','sub_charest','sub_esteban','sub_fabio','sub_gab','sub_gabriella','sub_genevieve','sub_gina','sub_jeremie','sub_julie','sub_katrine','sub_lisa'...
        ,'sub_marc','sub_marie','sub_mathieu','sub_maxime','sub_mingham','sub_patricia','sub_po','sub_russell','sub_sunachakan','sub_vincent'} ;
glm = load('c:/shared/gs') ; glm = mean(glm.gs) ; 
hrf = spm_hrf(2) ; %[0,1,6,9,7,4,-1,0] ; 
ntrs = 245 ; 
for sub=1:length(subs)
    cd(['c:/shared/allfmris/',subs{sub},'/ica']) ;  
    ics = load('melodic_mix') ; 
    braincomps = load_nii('melodic_IC.nii.gz') ; 
    bimg = braincomps.img ; 
    f1 = load_nii('f1.nii.gz') ; 
    f1im = f1.img ; 
    
    cd ../trigs ; clear alltrigs allparams
    alltimes = dir('stimTimes*.mat') ; 
    for stimFile=1:size(alltimes,1) ; % for all scans (1-9 usually)
        ctime = load(alltimes(stimFile).name) ;
        scanTimes = ctime.stimTimes ; 
        for presentation=1:size(scanTimes,2) % for all presentations within the scan
            alltrigs(stimFile,presentation) = scanTimes{presentation}(1) ; 
            allparams(stimFile,presentation) = scanTimes{presentation}(2) ; 
        end
    end
    stimes = round(squeeze(alltrigs(1,:))) ; % trigger presentation times
    trtimes = stimes/2 ; % converted to TR
    ideal = zeros(1,ntrs) ; 
    for i=1:size(trtimes,2) % fill the ideal with 1s where there was a stimulus
        ideal(trtimes(i)) = 1 ; 
    end
    idealhrf = conv(ideal,hrf) ; % convolve ideal with HRF
    idealhrf = idealhrf(1:245) ; % shave off last few indices
    glm = repmat(idealhrf,[1,9]) ; 
    
    g = graythresh(f1im) ; 
    mask = (mat2gray(f1im)>g/2) ; 
    %binmaskrows = mat2gray(sum(sum(mask,1),3)>0) ; maskinds = find(binmaskrows>0) ; binmaskrows(maskinds(1:10)) = 10 ; 
    %binmaskrows = smooth(binmaskrows) ; 

   % [specs,freqs] = spectopo(ics',0,0.5,'plot','off') ; 
   % stimfreqs = find(freqs>.06 & freqs<.065) ;
   % surroundfreqs = find((freqs>.07 & freqs <.15) | (freqs<.055 & freqs>.025)) ; 
    sumxz = squeeze(sum(sum(bimg(:,:,:,:),1),3)) ; 

    midpoint = round(size(bimg,2)/2) ; 
    for i=1:size(ics,2) ; stimhz(i) = corr2(ics(:,i)',glm) ; stimhz2(i) = corr2(ics(:,i)',glm); end
   % for i=1:size(ics,2) ; corrics(i) = corr2(binmaskrows,sumxz(:,i)) ; end
    %stimhz = mean(specs(:,stimfreqs),2) ; %-mean(specs(:,surroundfreqs),2) ; 
    for i=1:size(bimg,4) ; 
       cimg = imfilter(squeeze(bimg(:,:,:,i)),fspecial('gaussian',3,3)) ; 
       if max(max(max(cimg(:,1:midpoint,:)))) < max(max(max(cimg(:,midpoint:end,:))))
           stimhz(i) = -1 ;
       end
    end
    %[sv,si] = sort(stimhz.*corrics','descend') ;
    [sv,si] = sort(stimhz,'descend') ; 
    sortbimg = bimg(:,:,:,si) ; 
    % display the results of automatic IC selection
    zbars = squeeze(sum(sum(mask,1),2)) ; maxind = find(zbars==max(zbars)) ; 
    f = figure ; 
    for i=1:10 ; subplot(10,2,i) ; plot(squeeze(ics(1:240,si(i))),'LineWidth',2) ; title(['corr=',num2str(stimhz2(si(i))),'comp#=',num2str(i)]); end
    for i=1:10 ; subplot(3,10,20+i) ; imagesc(rot90(squeeze(mean(sortbimg(:,:,maxind-4:maxind,i),3)))) ; title(i) ;  end
    set(f, 'Position', [0, 0 1800 1000])
    suptitle(subs{sub}) ; 
    
    % save stuff for later analysis:
    subtrigs(sub,:,:) = alltrigs ; 
    subparams(sub,:,:) = allparams ; 
    submaps{sub} = squeeze(sortbimg(:,:,:,1)) ;
    subics(sub,:) = ics(:,si(1)) ; 
    
end

cd c:/shared/fmricomps
save('subtrigs','subtrigs') ; save('subparams','subparams') ; save('submaps','submaps') ; save('subics','subics') ; 





