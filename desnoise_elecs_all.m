clear all ; close all ; 
cd c:/shared/lastute/ ;
mongs=dir('*') ; mongs(1:2) = [] ; 
padamt = 15 ; 

for mong=1:length(mongs) ; 
    
cd(['C:\shared\lastute\',mongs(mong).name]) ; ls ; 

mask = load_untouch_nii('fnirt/ute_mask.nii.gz') ; 
resute = load_untouch_nii('res_ute.nii.gz') ; 
l = load_untouch_nii('layers.nii.gz') ; 

maskim = pad3d(mask.img > 0,padamt) ;
layerim = pad3d(l.img,padamt) ; 
layerim(layerim==0) = nan ; 
layerim = max(max(max(layerim)))-layerim ; 
layerim(isnan(layerim)) = 0 ; 
maskshell = imdilate(maskim,strel(ones(15,15,15))) - imerode(maskim,strel(ones(9,9,9))) ; 
uteim = pad3d(resute.img,padamt) ; 
uteim = uteim - imfilter(uteim,fspecial('gaussian',45,35)) ; 
smoothmask = imfilter(double(maskim),fspecial('gaussian',15,5)) ; 

shell2 = imdilate(maskim,strel(ones(11,11,11))) - imdilate(maskim,strel(ones(3,3,3))) ; 
elecint = shell2.*pad3d(resute.img,padamt) ; 
shellute = uteim.*double(maskshell) ; 
esize = 12 ; 

for eloc=1:20
    locs = load(['mricoords_',num2str(eloc),'.mat']) ; locs = locs.mricoords + padamt ; 
    clear finalint rawute
    for e=1:length(locs)
    loc1 = locs(:,e) ; 
    finalint(e,:,:,:) = squeeze(elecint(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize)) ; 
    rawute(e,:,:,:) = squeeze(shellute(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize)) ;
    end

    % do the denosing
    [gx,gy,gz] = ndgrid(-25:25,-25:25,-25:25) ; 
    [th,phi,r] = cart2sph(gx,gy,gz) ; 
    shell = r>6 & r<7 ; 
    inds = find(shell==1) ; 
    vecs  = [gx(inds),gy(inds),gz(inds)] ; 
    vecs = vecs./repmat(sum(vecs.^2,2),[1,3]) ; 
    blob = squeeze(finalint(1,:,:,:)) ; 
    blob_center = (size(squeeze(finalint(1,:,:,:))) + 1) / 2 ; 

    clear allinvs
    for ec=1:65 ; % figure,
        vecblobs = zeros(size(vecs,1),size(finalint,2),size(finalint,3),size(finalint,4)) ; 
        for v=1:length(vecs) ; 
            shellblob = squeeze(finalint(ec,:,:,:)) ; 
            rot = vrrotvec2mat(vrrotvec(vecs(v,:),[1,0,0])) ;
            rot(4,4) = 1 ; 
            T3 = [1 0 0 0
                  0 1 0 0
                  0 0 1 0
                  blob_center 1] ; 
            T1 = [1 0 0 0
                  0 1 0 0
                  0 0 1 0
                 -blob_center 1] ; 
            T = T1 * rot * T3 ; 
            tform = maketform('affine', T) ;
            tformfwd(blob_center, tform) ; 
            R = makeresampler('linear', 'fill') ; 
            TDIMS_A = [1 2 3];
            TDIMS_B = [1 2 3];
            TSIZE_B = size(blob);
            TMAP_B = [];
            F = 0;
            blob2 = tformarray(shellblob, tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F);
            vecblobs(v,:,:,:) = blob2 ; v,
        end
        maxvecblobs = max(vecblobs,[],4) ; 
        %icount =1 ; for i=1:10:380 ; subplot(6,7,icount) ; imagesc(squeeze(maxvecblobs(i,:,:))) ; icount = icount + 1 ; end 

        for i=1:size(maxvecblobs,1) ; repblobs(i,:,:,:) = repmat(squeeze(maxvecblobs(i,:,:)),[1,1,25]) ; end
        invblobs = zeros(size(vecs,1),size(finalint,2),size(finalint,3),size(finalint,4)) ; 
        for i=1:size(vecs,1)
            shellblob = squeeze(repblobs(i,:,:,:)) ; 
            rot = vrrotvec2mat(vrrotvec(vecs(i,:),[1,0,0])) ;
            rot(4,4) = 1 ; 
            rot = inv(rot) ; 
            T3 = [1 0 0 0
                0 1 0 0
                0 0 1 0
                blob_center 1] ; 
            T = T1 * rot * T3 ; 
            tform = maketform('affine', T) ;
            tformfwd(blob_center, tform) ; 
            R = makeresampler('linear', 'fill') ; 
            TDIMS_A = [1 2 3];
            TDIMS_B = [1 2 3];
            TSIZE_B = size(blob);
            TMAP_B = [];
            F = 0;
            blob2 = tformarray(shellblob, tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F);
            invblobs(i,:,:,:) = blob2 ; i,
        end
        allinvs(ec,:,:,:) = squeeze(mean(invblobs,1)) ; 
    end

    %for i=1:65 ; figure
    %    for j=1:17 ; subplot(5,5,j) ; imagesc(squeeze(allinvs(i,:,:,j)),[0,800]) ; end
    %end

    clear eleck eleck2 dilints
    for i=1:65 ; 
        ki = (squeeze(allinvs(i,:,:,:))) ; 
        [tv,topinds] = sort(ki(:),'descend') ; 
        newki = zeros(size(ki)) ; 
        newki(topinds(1:200)) = 1 ; 
        newki2 = zeros(size(ki)) ; 
        newki2(topinds(101:200)) = 1 ; 
        eleck(i,:,:,:) = newki ;
        eleck2(i,:,:,:) = newki2 ;   
        dilint = imdilate(newki,strel(ones(3,3,3))).*squeeze(rawute(i,:,:,:)) ; 
        inds = find(dilint~=0) ; [kc,km] = kmeanscustom(uint8(mat2gray(dilint(inds))*255),2) ; 
        dilint(inds) = km ; dilint = dilint==2 ; 
        dilints(i,:,:,:) = dilint ; 
    end
    %for i=1:65 ; figure ; for j=1:17 ; subplot(4,5,j) ; imagesc(squeeze(dilints(i,:,:,j))) ; end ; end

    segimg = zeros(size(shellute)) ; 
    boximg = zeros(size(shellute)) ; 
    colorimg = zeros(size(shellute)) ; 
    origimg = zeros(size(shellute)) ; 
    rawimg = zeros(size(shellute)) ; 
    
    for e=1:length(locs)
        loc1 = locs(:,e) ; 
        ek = squeeze(eleck(e,:,:,:)) ; 
        segimg(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize) = ek ; 
        boximg(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize) = 1 ; 
        colorimg(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize) = double(ek)*e ;           
        origimg(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize) = squeeze(allinvs(e,:,:,:)) ;               
        rawimg(loc1(1)-esize:loc1(1)+esize,loc1(2)-esize:loc1(2)+esize,loc1(3)-esize:loc1(3)+esize) = squeeze(rawute(e,:,:,:)) ;               
    end
    
    resute.img = segimg(padamt:end-(padamt+1),padamt:end-(padamt+1),padamt:end-(padamt+1)) ; 
    resute.img = uint8(resute.img) ; 
    save_untouch_nii(resute,['resute_seg',num2str(eloc),'.nii.gz']) ; 

    resute.img = boximg(padamt:end-(padamt+1),padamt:end-(padamt+1),padamt:end-(padamt+1)) ; 
    save_untouch_nii(resute,['boximg',num2str(eloc),'.nii.gz']) ; 

    resute.img = colorimg(padamt:end-(padamt+1),padamt:end-(padamt+1),padamt:end-(padamt+1)) ; 
    save_untouch_nii(resute,['colorimg',num2str(eloc),'.nii.gz']) ; 
    
    resute.img = origimg(padamt:end-(padamt+1),padamt:end-(padamt+1),padamt:end-(padamt+1)) ; 
    save_untouch_nii(resute,['origimg',num2str(eloc),'.nii.gz']) ; 
    
    resute.img = rawimg(padamt:end-(padamt+1),padamt:end-(padamt+1),padamt:end-(padamt+1)) ; 
    save_untouch_nii(resute,['rawimg',num2str(eloc),'.nii.gz']) ; 
end
end




