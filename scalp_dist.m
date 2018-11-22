clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 

for sb=1:length(subs); disp(sb) 
    cd(['c:/shared/freesurfer_segs/sub_',subs{sb},'/mri']) ;  
    ctx = load_untouch_nii('cortex.nii.gz');
    
    %[cort_x,cort_y,cort_z] = ind2sub(size(ctx.img),find(ctx.img>0)); 
    
    ctximg = double(imdilate(ctx.img,strel(ones(5,5,5)))); 

    cd(['c:/shared/allfmris/sub_',subs{sb}]) ;  
    brain = load_untouch_nii('fs_brain.nii.gz'); 
    tcorrs = load_untouch_nii('cleancorrs_fs.nii.gz'); 
    tcorrs.img(:,:,100:end) = 0; 
    %tcorrs.img = imfilter(double(tcorrs.img),fspecial('gaussian',11,9)); 
    [sv,si] = sort(tcorrs.img(:),'descend'); 
    zcorrs = zeros(size(tcorrs.img));
    zcorrs(si(1:10000)) = 1; 
    zcorrs = double(zcorrs).*tcorrs.img; 
    
    %[cort_x,cort_y,cort_z] = ind2sub(size(zcorrs),find(zcorrs>0)); 
    
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
        
    
    
end
cd(['E:\clean_allres\vincent']) ; ls 
alltopos = load('alltopos'); alltopos = alltopos.alltopos; 

cd C:\shared\allres\alex
eeg = pop_loadset('resamp_vis09.set'); 
elabs = {eeg.chanlocs.labels}; 
cd c:/shared/ ; 
elecorder = load('elecorder.mat'); elecorder = elecorder.elecorder; 
for i=1:length(elecorder)
    eleci = elecorder{i};
    indi = find(strcmpi(eleci,elabs)); 
    if ~isempty(indi)
        sbothdists(:,indi) = comdists(:,i); 
        
    end
end

for i=1:64
    subplot(5,13,i);
   plot(squeeze(sbothdists(:,i)),squeeze(alltopos(:,i)),'o') ;
    hl = lsline; 
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    Slope = B(2); 
    allslopes(i) = Slope; 
    title([num2str(Slope)]);
end

for i=1:64 ; [ecorrs(i),eps(i)] = corr(sbothdists(:,i),alltopos(:,i)) ; end
figure,
subplot(1,2,2) ; 
[sv,si] = sort(ecorrs,'ascend'); 
for i=1:15 ; plot(sbothdists(:,si(i)),alltopos(:,si(i)),'d','LineWidth',6) ;  hold on ; end ;% ylim([-0.01,0.12]) ; xlim([20,140]); 
legend(elabs(si(1:15)))
subplot(1,2,1) ; topoplot(eps,eeg.chanlocs,'maplimits',[0,0.1],'electrodes','labels'); colormap cool ; title('R-squared')


figure,plot(mean(alltopos,1),ecorrs.^2,'d','LineWidth',6); lsline; 
[r,p] = corr(mean(alltopos(:,~isnan(ecorrs)),1)',(ecorrs(~isnan(ecorrs)).^2)');
title(['r=',num2str(r)]); 

badchans = [28,32,41,17,46,22]; 
goodchans = zeros(1,size(sbothdists,2));
goodchans(badchans) = 1; goodchans = find(goodchans==0); 
postelecs = [23,56,24,57,25,58,26,59,27,60,61,62,63,64,29,30,31];
[gp_sv,gp_si] = sort(mean(alltopos(:,postelecs),2),'ascend'); 
colors = colormap(parula); colors = imresize(colors,[24,3]); colors(colors>1) = 1; colors(colors<0) = 0; 
[subsort_val,subsort_ind] = sort(mean(alltopos(:,postelecs),2),'descend'); 


subplot(2,8,5);
topoplot(allslopes,eeg.chanlocs) ; title('slope'); colorbar; 
subplot(2,8,6) ;plot(mean(sbothdists(:,goodchans),1),abs(allslopes(goodchans)),'ko','LineWidth',1) ; lsline; xlabel('electrode to ROI distance (mm)'); ylabel('slope'); 
[r,p] = corr(mean(sbothdists(:,goodchans),1)',allslopes(goodchans)'); 
title(['r=',num2str(r),', p=',num2str(p)]); 

subplot(2,8,1); 
plot(mean(sbothdists(:,postelecs),2),mean(alltopos(:,postelecs),2),'ko','LineWidth',1) ; lsline; xlim([50,64]); ylim([-0.01,.1]); xlabel('distance');  ylabel('gamma amplitude'); 
[c,p] = corr(mean(sbothdists(:,postelecs),2),mean(alltopos(:,postelecs),2)); 
title(['rho=',num2str(c),', p=',num2str(p)]);

for i=1:64 ; [chancorrs(i),chanps(i)]  =corr(sbothdists(:,i),alltopos(:,i)) ; end
subplot(2,8,2),topoplot(chancorrs,eeg.chanlocs,'style','map'); title('rho') ; colorbar; colormap parula; 
subplot(2,8,3); topoplot(chanps,eeg.chanlocs,'maplimits',[0,0.1],'style','map'); colorbar ; title('p');colormap parula;
postlabs = elabs(postelecs); 
subplot(2,8,4) ; barh(chancorrs(postelecs)); set(gca,'YTick',1:length(postelecs),'YTickLabel',postlabs); title('correlation at all posterior electrodes'); vline(-0.4,'r'); text(-0.4,19,'p<0.05 (n-24)'); 
xlabel('correlation'); ylabel('electrode'); 

[sv,sort_elecinds] = sort(mean(sbothdists(:,postelecs),2),'descend'); 
subplot(2,8,7) ; bar(mean(alltopos(sort_elecinds,postelecs),2)); title('sorted gamma amplitude'); xlabel('high dist <......> low dist') ; ylabel('gamma amplitude'); 
[h,p,ci,stats] = ttest2(mean(alltopos(sort_elecinds(1:10),postelecs),2),mean(alltopos(sort_elecinds(end-9:end),postelecs),2)); 
topsubs = mean(squeeze(alltopos(sort_elecinds(1:10),postelecs)),2); botsubs = mean(squeeze(alltopos(sort_elecinds(end-9:end),postelecs)),2); diffsubs = [topsubs,botsubs]; 
tdist = mean(mean(sbothdists(sort_elecinds(1:10),postelecs),1),2); bdist = mean(mean(sbothdists(sort_elecinds(end-10:end),postelecs),1),2); 

subplot(2,8,8) ; 
barwitherr(std(diffsubs,0,1)/sqrt(10),mean(diffsubs,1)); xlim([0.5,2.5]); ylim([0,.065]); 
perc = round(((mean(diffsubs(:,2))-mean(diffsubs(:,1)))/mean(diffsubs(:,1)))*100); 
title(['p=', num2str(p),', percent increase = ',num2str(perc),'%']); 
set(gca,'XTickLabel',{'high distance','low distance'}); ylabel('gamma amplitude'); xlabel('group (mean difference = 5mm)'); 


mydata=mean(alltopos(sort_elecinds,postelecs),2);
sortdata = mean(alltopos(gp_si,postelecs),2); 
subplot(2,8,10)
hold on
for i=1:length(mydata)
    h=bar(i,mydata(i));
    ind = find(sort(mydata)==mydata(i));     
    set(h,'FaceColor',[colors(ind,1),colors(ind,2),colors(ind,3)]);    
end
set(gca,'XTick',1:5:24,'XTickLabel',round(sv(1:5:24))); xlabel('mean distance (mm)'); ylabel('gamma amplitude'); xlim([0.5,24.5]); 
title('amplitude sorted by distance'); 
subplot(2,8,9); hold on
for i=1:length(mydata)
    h=bar(i,sortdata(((i))));
    set(h,'FaceColor',[colors(i,1),colors(i,2),colors(i,3)]);    
end
title('amplitude sorted by amplitude'); xlabel('subject') ; ylabel('gamma amplitude'); xlim([0.5,24.5]); 

subplot(2,8,11) ; 
h = bar(1,mean(diffsubs(:,1),1)); hold on ;
set(h,'FaceColor',[colors(1,1),colors(1,2),colors(1,3)]);    
h = bar(2,mean(diffsubs(:,2),1)); 
set(h,'FaceColor',[colors(24,1),colors(24,2),colors(24,3)]);    
errorbar(mean(diffsubs,1),std(diffsubs,0,1)/sqrt(10),'k.'); 
xlim([0.5,2.5]); ylim([0,.065]); 
title(['p=', num2str(p),', percent increase = ',num2str(perc),'%']); 
set(gca,'XTick',[1,2],'XTickLabel',{['mean dist=',num2str(round(tdist)),'mm'],['mean dist=',num2str(round(bdist)),'mm']}); ylabel('gamma amplitude'); 

subplot(2,8,12) ; topoplot(mean(alltopos,1),eeg.chanlocs) ; colormap parula ; colorbar;
subplot(2,8,13) ; imagesc([min(mean(alltopos(:,postelecs),2)),max(mean(alltopos(:,postelecs),2))]); colorbar; 

% plot all electrodes, scatter plot according to distance
subplot(2,8,14); 
newcolors = imresize(colors,[17,3]); newcolors(newcolors>1) = 1; newcolors(newcolors<0) = 0; 
for i=1:length(postelecs)
   plot(squeeze(sbothdists(:,postelecs(i))),squeeze(alltopos(:,postelecs(i))),'d','LineWidth',3,'Color',[newcolors(18-i,1),newcolors(18-i,2),newcolors(18-i,3)]); hold on ;   
end
xlabel('distance(mm)'); ylabel('gamma amplitude'); 
subplot(2,8,15);
imagesc([min(mean(sbothdists(:,postelecs))),max(mean(sbothdists(:,postelecs)))]); colorbar ; 




postchans = [60,61,62,63,64,29,30,31,23,56,24,57,25,58,26,59,27];

%mtmersp = squeeze(mean(mean(allmersp(:,:,:,40:160),2),4)); 

[corrs,ps] = corr(sbothdists,alltopos); 
[mcorrs,mps] = corr(alltopos,mean(sbothdists(:,postchans),2)); 


for i=1:64 ; subplot(5,13,i);
    topoplot((corrs(:,i)),eeg.chanlocs,'maplimits',[-.8,.8],'plotchans',goodchans)

end

figure,
subplot(2,2,1); 
plot(mcorrs) ; ylim([-.6,.6]) ; hline([0.4,-0.4],'k'); hline(0,'r');  



figure,
subplot(1,2,1);
topoplot(squeeze(mean(corrs(:,20:40),2)),eeg.chanlocs,'plotchans',goodchans,'electrodes','labels','maplimits',[-.55,.55])
subplot(1,2,2); 
topoplot(squeeze(mean(sbothdists(:,:),1)),eeg.chanlocs,'plotchans',goodchans,'electrodes','labels','maplimits',[0,150])









