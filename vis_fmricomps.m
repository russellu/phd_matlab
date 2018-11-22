clear all ; close all ; 
subs = {'sub_alex','sub_charest','sub_esteban','sub_fabio','sub_gab','sub_gabriella','sub_genevieve','sub_gina','sub_guillaume','sub_jeremie','sub_julie',...
'sub_katrine','sub_lisa','sub_marc','sub_marie','sub_mathieu','sub_maxime','sub_mingham','sub_patricia','sub_po','sub_russell','sub_sunachakan','sub_tah','sub_vincent'} ;
for sb=1:length(subs)
cd(['c:/shared/allfmris/',subs{sb}]) ;
gm = load_untouch_nii('epi_dilcortex.nii.gz') ; 
commons = dir('mel_common*') ; 
for i=1:length(commons)
    cd(['c:/shared/allfmris/',subs{sb},'/',commons(i).name]) ; 
    meanfmri = load_untouch_nii('mean.nii.gz') ;
    ics = load_untouch_nii('melodic_IC.nii.gz') ; 
    mix = load('melodic_mix') ; 
    f = abs(fft(mix,[],1)) ; alff = sum(f(2:14,:)) - sum(f(20:50,:),1) ; 
    %figure,
    for j=1:size(ics.img,4) ; 
    %    subplot(7,10,j) ; 
        imgj = ics.img(:,:,:,j) ; 
        [sv,si] = sort(imgj(:),'descend') ; 
        zimg = zeros(size(imgj)) ; zimg(si(1:2000)) = 1 ; 
        ximg = zimg.*gm.img ; perc = sum(ximg(:))/2000 ; 
     %   imagesc(squeeze(max(ics.img(:,:,:,j),[],3))) ; colormap jet ; title(['perc= ',num2str(perc),' n=',num2str(j)]) ; 
      %  set(gca,'XTick',[],'YTick',[]) ; 
        percs(j) = perc ;      
    end
    goodcs = find(percs>.6) ; 
    badcs = zeros(1,size(mix,2)) ; badcs(goodcs) = 1 ; bads = find(badcs==0) ; 
    dlmwrite('badcs.txt',bads) ; 
end
end

