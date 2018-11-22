clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 
elabs = load('C:\shared\all_white_normals\a2_good\elabs') ; elabs = elabs.elabs ; 
elecorder = load('C:\shared\all_white_normals\a2_good\elecorder') ; elecorder = elecorder.elecorder ; 
baseeg = pop_loadset('c:/shared/allres/alex/cleanfilt.set') ; 
baselabs = {baseeg.chanlocs.labels} ;
for i=1:length(baselabs) ;
    if ~isempty(find(strcmpi(baselabs{i},elecorder)))
        orderinds(i) = find(strcmpi(baselabs{i},elecorder)) ;
    end
end 

for sb=1:length(subs); disp(sb) ; 
    cd(['c:/shared/allfmris/sub_',subs{sb}]) ;  

    locs = load('locs') ; locs = locs.locs ; 
    meancorrs = load_untouch_nii('trim_cleancorrs_fs.nii.gz') ; t1 = load_untouch_nii('fs_t1.nii.gz') ; 
    
    threshs = 0.4:0.03:0.99 ; 
    for thr=1:length(threshs)
    threshcorrs = meancorrs.img>threshs(thr) ; sumcorrs(sb,thr) = sum(threshcorrs(:)) ; 
    end
    [sv,si] = sort(meancorrs.img(:),'descend') ; 
    zcorrs = zeros(size(meancorrs.img)) ; zcorrs(si(1:8000)) = 1 ; 
    [sx,sy,sz] = ind2sub(size(meancorrs.img),si(1:8000)) ; 
    [cx,cy,cz] = centmass3(zcorrs) ; 
    sqrdiffs = sqrt((locs(1,:)-cx).^2 + (locs(2,:)-cy).^2 + (locs(3,:)-cz).^2) ; 
    for i=1:size(locs,2)
        for j=1:length(sx)
            mindiffs(i,j) = sqrt((locs(1,i)-sx(j)).^2 + (locs(2,i)-sy(j)).^2 + (locs(3,i)-sz(j)).^2) ;
        end
    end
    mindiffs = min(mindiffs,[],2) ;
    
    distmat = zeros(1,64) ; 
    for i=1:length(orderinds) ; if orderinds(i)~=0 ; distmat(i) = sqrdiffs(orderinds(i)) ; end ; end
    for i=1:length(orderinds) ; if orderinds(i)~=0 ; mindistmat(i) = mindiffs(orderinds(i)) ; end ; end

    bads = find(distmat==0) ; goods = find(distmat~=0) ; 
    goods(find(goods==17)) = [] ; goods(find(goods==22)) = [] ; goods(find(goods==41)) = [] ; goods(find(goods==46)) = [] ; 
    alldistmat(sb,:) = distmat ; allmindistmat(sb,:) = mindistmat ; 
    cd(['c:/shared/allres/',subs{sb}]) ; 
    alldiffs = load('alldiffs.mat') ; alldiffs = alldiffs.alldiffs ; alleeg(sb,:,:,:) = alldiffs ; 

end
for i=1:size(alleeg,1) ; for j=1:size(alleeg,2) ; for k=1:size(alleeg,3) ; smootheeg(i,j,k,:) = imfilter(squeeze(alleeg(i,j,k,:)),fspecial('gaussian',[9,1],5)) ; end ; end ; end
postelecs = [23,56,24,57,25,58,26,59,27,60,61,62,63,64,29,30,31] ;


