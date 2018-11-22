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
for sb=1:length(subs) ; disp(sb) ; 
    cd(['c:/shared/allfmris/sub_',subs{sb},'']) ;  
    locs = load('locs') ; locs = locs.locs ; 
    meancorrs = load_untouch_nii('cleancorrs_fs.nii.gz') ; 
    meancorrs.img(:,:,110:end) = 0 ; 
    [sv,si] = sort(meancorrs.img(:),'descend') ; 
    zcorrs = zeros(size(meancorrs.img)) ; zcorrs(si(1:5000)) = 1 ; 
    threshs = 0.4:0.03:0.99 ; 
    for thr=1:length(threshs)
        threshcorrs = meancorrs.img>threshs(thr) ; 
        sumcorrs(sb,thr) = sum(threshcorrs(:)) ; 
    end
    cd(['c:/shared/gamma_t1s/sub_',subs{sb},'/mri']) ;  
    leftnorms = load_untouch_nii('surf.nii.gz') ; 
    rightnorms = load_untouch_nii('surfr.nii.gz') ; 
    bothnorms = leftnorms.img + rightnorms.img ; 
    
    normls = 100:500:15000 ; 
    for nl=1:length(normls) ; 
    [nx,ny,nz] = ind2sub(size(zcorrs),si(1:normls(nl))) ; 
    allnorms = zeros(length(nx),3) ; 
    for i=1:length(nx)
        allnorms(i,:) = squeeze(bothnorms(nx(i),ny(i),nz(i),:)) ; 
    end
    inots(sb,nl) = 1 - sum(abs(sum(allnorms,1)))./length(nx) ; 
    end
    cd(['c:/shared/allres/',subs{sb}]) ; 
    alldiffs = load('alldiffs.mat') ; alldiffs = alldiffs.alldiffs ; alleeg(sb,:,:,:) = alldiffs ;   

end

for i=1:size(alleeg,1) ; for j=1:size(alleeg,2) ; for k=1:size(alleeg,3) ; smootheeg(i,j,k,:) = imfilter(squeeze(alleeg(i,j,k,:)),fspecial('gaussian',[9,1],5)) ; end ; end ; end
meeg = squeeze(mean(smootheeg(:,:,:,:),2)) ; 
postelecs = [23,56,24,57,25,58,26,59,27,60,61,62,63,64,29,30,31] ;
emeeg=  squeeze(mean(meeg(:,postelecs,:),2)) ; 
[c,p] = corr(emeeg,inots)  ;
[cs,ps] = corr(emeeg,sumcorrs) ; 
reg1 = inots(:,10) ; reg2 = sumcorrs(:,10) ; 
regs = [reg1,reg2] ; 
mmeeg = squeeze(mean(meeg(:,postelecs,:),2)) ; 
for i=1:128
b(:,i) = mvregress(regs,mmeeg(:,i)) ; 
end
for i=1:128 
    corrs(i) = corr2(mmeeg(:,i),b(1,i)*reg1+b(2,i)*reg2) ; 
end






