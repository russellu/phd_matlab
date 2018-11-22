clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 

for sb=1:length(subs); disp(sb) 
    cd(['c:/shared/freesurfer_segs/sub_',subs{sb},'/mri']) ;  
    ctx = load_untouch_nii('cortex.nii.gz');    
    ctximg = double(imdilate(ctx.img,strel(ones(5,5,5)))); 

    cd(['c:/shared/allfmris/sub_',subs{sb}]) ;  
    brain = load_untouch_nii('fs_brain.nii.gz'); 
    tcorrs = load_untouch_nii('cleancorrs_fs.nii.gz'); 
    tcorrs.img(:,:,100:end) = 0; 
    [sv,si] = sort(tcorrs.img(:),'descend'); 
    zcorrs = zeros(size(tcorrs.img));
    zcorrs(si(1:10000)) = 1; 
    zcorrs = double(zcorrs).*tcorrs.img; 
        
    [cort_x,cort_y,cort_z] =   ind2sub(size(zcorrs),find(zcorrs>0));
    [cort_x2,cort_y2,cort_z2] =  centmass3((zcorrs)); 

       
    locs = load('locs.mat'); 
    locs = locs.locs; 
    
    for i=1:size(locs,2)
        dists_i = sqrt((cort_x-locs(1,i)).^2 + (cort_y-locs(2,i)).^2 + (cort_z-locs(3,i)).^2); 
        mindists(sb,i) = min(dists_i); 
        comdists(sb,i) = sqrt((cort_x2-locs(1,i)).^2 + (cort_y2-locs(2,i)).^2 + (cort_z2-locs(3,i)).^2); 
    end
    cd(['c:/shared/allres/',subs{sb}]) ;  
    mersp = load('amersp') ; mersp = mersp.amersp ; allmersp(sb,:,:,:) = mersp ;     
    
    
    %inot
    bininds = find((tcorrs.img(:)>.6)); 
    cd(['c:/shared/gamma_t1s/sub_',subs{sb},'/mri']) ;  
    leftnorms = load_untouch_nii('surf.nii.gz') ; 
    rightnorms = load_untouch_nii('surfr.nii.gz') ; 
    bothnorms = leftnorms.img + rightnorms.img ; 
    normls = 100:500:25000 ; 
    for nl=1:length(normls)  
    [nx,ny,nz] = ind2sub(size(zcorrs),si(1:normls(nl))) ; 
    [bx,by,bz] = ind2sub(size(zcorrs),bininds); 
    allnorms = zeros(length(nx),3) ; 
    allbnorms = zeros(length(bx),3); 
    for i=1:length(nx)
        allnorms(i,:) = squeeze(bothnorms(nx(i),ny(i),nz(i),:)) ; 
    end
    for i=1:length(bx)
        allbnorms(i,:) = squeeze(bothnorms(bx(i),by(i),bz(i),:)) ; 
    end
    binots(sb) = 1 - sum(abs(sum(allbnorms,1)))./length(bx) ; 
    inots(sb,nl) = 1 - sum(abs(sum(allnorms,1)))./length(nx) ; 
    end
end
cd(['E:\clean_allres\vincent']) ; ls 
alltopos = load('alltopos'); alltopos = alltopos.alltopos; 

postelecs = [23,56,24,57,25,58,26,59,27,60,61,62,63,64,29,30,31];
[postcorrs,postps] = corr(squeeze(mean(alltopos(:,postelecs),2)),inots); 
[meancorr,meanp] = corr(squeeze(mean(alltopos(:,postelecs),2)),mean(inots,2)); 
[eleccorrs,elecps] = corr(alltopos,inots); 
[sv_inot,si_inot] = sort(mean(inots,2),'descend'); 
m_inots = mean(inots,2); 
sort_inots = [mean(alltopos(si_inot(1:10),postelecs),2),mean(alltopos(si_inot(end-9:end),postelecs),2)]; 

subplot(2,8,1); 
plot(squeeze(mean(inots,2)),squeeze(mean(alltopos(:,postelecs),2)),'ok'); lsline;
xlabel('I0') ; ylabel('gamma amplitude'); 
[c,p] = corr(squeeze(mean(inots,2)),squeeze(mean(alltopos(:,postelecs),2)));
title(['rho=',num2str(c),' p=',num2str(p)]);

subplot(2,8,3); 
colors = colormap(parula); colors = imresize(colors,[24,3]); colors(colors>1) = 1; colors(colors<0) = 0; 
perc = round(((mean(sort_inots(:,2))-mean(sort_inots(:,1)))/mean(sort_inots(:,1)))*100); 
[~,p,ci,stats] = ttest2(sort_inots(:,1),sort_inots(:,2)); 
h = bar(1,mean(sort_inots(:,1),1)); hold on ;
set(h,'FaceColor',[colors(1,1),colors(1,2),colors(1,3)]);    
h = bar(2,mean(sort_inots(:,2),1)); 
set(h,'FaceColor',[colors(24,1),colors(24,2),colors(24,3)]);    
errorbar(mean(sort_inots,1),std(sort_inots,0,1)/sqrt(10),'k.'); 
xlim([0.5,2.5]); ylim([0,.065]); 
title(['p=', num2str(p),', percent increase = ',num2str(perc),'%']); 
set(gca,'XTick',[1,2],'XTickLabel',{['mean I0=',num2str(round(mean(m_inots(si_inot(1:10)))*100)/100)],['mean I0=',num2str(round(mean(m_inots(si_inot(end-9:end)))*100)/100)]}); ylabel('gamma amplitude'); 

subplot(2,8,2); 
mydata = squeeze(mean(alltopos(si_inot,postelecs),2)); 
hold on
for i=1:length(mydata)
    h=bar(i,mydata(i));
    ind = find(sort(mydata)==mydata(i));     
    set(h,'FaceColor',[colors(ind,1),colors(ind,2),colors(ind,3)]);    
end
set(gca,'XTick',1:5:24,'XTickLabel',round(sv_inot(1:5:24)*1000)/1000); xlabel('i0'); ylabel('gamma amplitude'); xlim([0.5,24.5]); 
title('amplitude sorted by I0'); 


