clear all ; close all ; 
subs= {'alex','dina','genevieve','jeremie','karl','russell','tegan','valerie'} ;

for s=1:length(subs) ; 
    clear ts tsinds corrinds maxcorrs
    cd(['c:/shared/badger_mri/',subs{s},'/nii']) ; ls 
    allangles = load_untouch_nii('allangles.nii.gz') ; 
    img = allangles.img ; 

    % get the average tuning curve
    % first get the voxels that responded to the stimulus (that are tuned)
    rsindex = max(img,[],4) - min(img,[],4) ; 

    resimg = reshape(rsindex,[1,numel(rsindex)]) ; 
    nz = find(resimg~=0) ; 
    meanthresh = mean(resimg(nz))./2 ; 
    tsinds = find((rsindex > meanthresh)>0) ; 
    [cx,cy,cz] = ind2sub(size(rsindex),find((rsindex > meanthresh)>0)) ; 
    for i=1:length(cx)
       ts(i,:) = squeeze(img(cx(i),cy(i),cz(i),:)) ; 
    end

    ts = ts./repmat(std(ts,0,2),[1,90]) ; 

    %{
    km = kmeans(ts(1:500),5) ; clear allt1s ; 
    for i=1:5
    k1 = find(km==i) ; 
    t1 = ts(k1,:) ; 
    subplot(4,5,i) ; imagesc(t1) ; 
    allt1s{i} = t1 ; 
    end
    refsig = mean(allt1s{3}) ; 
    %}
    xs = -45:44 ; 
    gx = exp(-(xs.^2)./500) ; refsig = circshift(gx',45) ; 

    clear corrs 
    for i=1:size(ts,1)
        jcount = 1 ; 
        for j=1:90
            newts = circshift(ts(i,:)',j) ; 
            corrs(i,jcount) = corr2(newts,refsig) ; 
            jcount = jcount + 1 ; 
        end
    end

    for i=1:size(corrs,1) ; corrinds(i) = find(corrs(i,:)==max(corrs(i,:)),1) ; maxcorrs(i) = max(corrs(i,:)) ; end 
    ref = load_untouch_nii('fref.nii.gz') ; 
    newimg = zeros(size(rsindex)) ; 
    newimg(tsinds) = corrinds ; 
    ref.img = newimg ; save_untouch_nii(ref,'retindices.nii.gz') ; 
    corrimg = zeros(size(rsindex)) ; 
    corrimg(tsinds) = maxcorrs ; 
    ref.img = corrimg ; save_untouch_nii(ref,'retcorrs.nii.gz') ; 


end

