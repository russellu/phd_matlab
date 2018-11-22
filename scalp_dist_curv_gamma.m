clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 

badsubs = [7,10,21]; goodsubs = zeros(1,24); goodsubs(badsubs)  = 1; goodsubs = find(goodsubs==0); 
peaksubs = [2,3,4,5,6,7,8,10,11,12,13,15,17,18,20,21,22,23] ; 
allpeaks = [64,76,74,64,66,84,72,66,74,74,70,66,64,70,64,88,70,56,84,60,66,66,54,84] ; 
for sb=1:length(subs); disp(sb) 
    
    % DISTANCE
    cd(['c:/shared/allfmris/sub_',subs{sb}]) ;  
    rh_ribbon = load_untouch_nii('rh_ribbon.nii.gz');
    lh_ribbon = load_untouch_nii('lh_ribbon.nii.gz'); 
    brain = load_untouch_nii('fs_brain.nii.gz'); 
    tcorrs = load_untouch_nii('cleancorrs_fs.nii.gz'); 
    tcorrs.img(:,:,100:end) = 0; 
    [sv,si] = sort(tcorrs.img(:),'descend'); 
    zcorrs = zeros(size(tcorrs.img));
    zcorrs(si(1:10000)) = 1; 
    zcorrs = double(zcorrs).*tcorrs.img; 
    sum_tcorrs(sb) = sum(tcorrs.img(:)>0.5); 
    right_zcorrs = zcorrs.*double(rh_ribbon.img);
    left_zcorrs = zcorrs.*double(lh_ribbon.img);
    sum_rightcorrs(sb) = sum(right_zcorrs(:)>0);
    sum_leftcorrs(sb) = sum(left_zcorrs(:)>0); 
    
    [cort_x,cort_y,cort_z] =   ind2sub(size(zcorrs),find(zcorrs>0));
    [cort_x2,cort_y2,cort_z2] =  centmass3((zcorrs)); 
    
    [right_cort_x,right_cort_y,right_cort_z] =  centmass3((right_zcorrs)); 
    [left_cort_x,left_cort_y,left_cort_z] =  centmass3((left_zcorrs)); 
        
    locs = load('locs.mat'); 
    locs = locs.locs; 
    
    for i=1:size(locs,2)
        dists_i = sqrt((cort_x-locs(1,i)).^2 + (cort_y-locs(2,i)).^2 + (cort_z-locs(3,i)).^2); 
        mindists(sb,i) = min(dists_i); 
        comdists(sb,i) = sqrt((cort_x2-locs(1,i)).^2 + (cort_y2-locs(2,i)).^2 + (cort_z2-locs(3,i)).^2); 
        
        right_dists(sb,i) = sqrt((right_cort_x-locs(1,i)).^2 + (right_cort_y-locs(2,i)).^2 + (right_cort_z-locs(3,i)).^2); 
        left_dists(sb,i) = sqrt((left_cort_x-locs(1,i)).^2 + (left_cort_y-locs(2,i)).^2 + (left_cort_z-locs(3,i)).^2);         
    end
    
    % I-NOT
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
    
    % ERSP
    cd(['c:/shared/allres/',subs{sb}]) ;  
    %mersp = load('amersp') ; mersp = mersp.amersp ; 
    %allmersp(sb,:,:,:) = mersp ;                    
end

cd C:\shared\allfmris\sub_gabriella
tcorrs = load_untouch_nii('cleancorrs.nii.gz'); 
cd trigs ; 
stimes = dir('stimTimes*mat') ; 
for st=1:length(stimes) ; sti = load(stimes(st).name) ; sti = sti.stimTimes ; sti = cell2mat(sti) ; allstims(st,:) = sti(1:2:end) ; end 
meantimes = mean(allstims,1) ; 
vec = zeros(1,490) ; for i=1:length(meantimes) ; vec(round(meantimes(i)):round(meantimes(i))+1) = 1 ; end 
conved = conv(vec,spm_hrf(1),'full') ; conved = conved(1:length(vec)) ; 
conved = imresize(conved,[1,245]) ; vec2 = imresize(vec,[1,245],'nearest'); 
cd .. ;
denoise_fmri = load_untouch_nii('denoise.nii.gz'); 
tcorrinds = find(tcorrs.img>0.7); resdenoise = reshape(denoise_fmri.img,[66*67*45,245]); 

cd(['E:\clean_allres\vincent']) ; ls 
mampdiff = load('mampdiff'); mampdiff = mampdiff.mampdiff; 
stimampdiff = mampdiff;
mampdiff = squeeze(mean(mampdiff,2)); 
submersp = load('submersp'); allmersp = submersp.submersp; 
times = load('times'); times = times.times;
freqs = load('freqs'); freqs = freqs.freqs; 
postchans = [60,61,62,63,64,29,30,31,23,56,24,57,25,58,26,59,27];
chanampdiff = squeeze(mean(mampdiff(:,postchans,:),2)); 
mersp = squeeze(mean(mean(mean(allmersp(:,:,20:50,40:80),2),3),4)); 
imagesc([],1:2:120,squeeze(mean(mean(allmersp))),[-2,.5]) ; axis xy ; hline([95,35],'k');hline([6,24],'r'); 
badchans = [28,32]; goodchans = zeros(1,64); goodchans(badchans) = 1; goodchans = find(goodchans==0); 

allstds = load('allstds'); allstds = allstds.allstds; stds = allstds; 
substds = (squeeze(mean(mean(mean(stds(:,:,1:2,:,20:end,:),3),4),5)));
[sv,si] = sort(substds,3,'descend'); 
allampdiff = load('allampdiff'); allampdiff = allampdiff.allampdiff; ampdiff = allampdiff; 
smampdiff = zeros(size(ampdiff,1),size(ampdiff,2),size(ampdiff,3),size(ampdiff,4)); 
gammastds = zeros(24,6); 
alphastds = zeros(24,6); 
cleanamps = zeros(24,64,50,120*6); 
for i=1:24
    for j=1:6
        smampdiff(i,j,:,:) = squeeze(mean(ampdiff(i,j,:,:,si(i,j,10:end)),5)); 
        gammastds(i,j) = squeeze(std(mean(mean(ampdiff(i,j,postchans,20:50,si(i,j,10:end)),3),4),0,5)); 
        cleanamps(i,:,:,(j-1)*120+1:j*120) =  squeeze(ampdiff(i,j,:,:,si(i,j,11:end))); 
    end
end

% get the significant electrodes in each subject/frequency band
melec_alpha = squeeze(mean(mean(cleanamps(:,:,4:12,:),3),4)); 
melec_gamma = squeeze(mean(mean(cleanamps(:,:,20:40,:),3),4)); 
for i=1:64
    [h,p,ci,t] = ttest(melec_gamma(:,i));
    stats_gamma(i) = t.tstat; 
    [h,p,ci,t] = ttest(melec_alpha(:,i));
    stats_alpha(i) = t.tstat; 
end

gamma_topo = squeeze(mean(mean(mean(cleanamps(:,:,20:50,:),1),3),4));
alpha_topo = squeeze(mean(mean(mean(cleanamps(:,:,4:12,:),1),3),4));

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
        srightdists(:,indi) = right_dists(:,i);
        sleftdists(:,indi) = left_dists(:,i); 
    end
end
for i=1:64
    for j=1:50
        ecorrs(i,j) = corr(squeeze(mampdiff(:,i,j)),squeeze(sbothdists(:,i))); 
    end
end
fcorrs = corr(squeeze(mean(mampdiff(:,postchans,:),2)),mean(sbothdists(:,postchans),2)); 
for i=1:64
[gamma_ecorrs(i),gamma_eps(i)] = corr(sbothdists(:,i),squeeze(mean(mampdiff(:,i,20:50),3))); 
[r_gamma_ecorrs(i),r_gamma_eps(i)] = corr(srightdists(:,i),squeeze(mean(mampdiff(:,i,20:50),3))); 
[l_gamma_ecorrs(i),l_gamma_eps(i)] = corr(sleftdists(:,i),squeeze(mean(mampdiff(:,i,20:50),3))); 
[alpha_ecorrs(i),alpha_eps(i)] = corr(sbothdists(:,i),squeeze(mean(mampdiff(:,i,4:12),3))); 
end

mean_amp = squeeze(mean(mean(mampdiff(:,postchans,20:50),2),3)); 
[sv_mean_amp,si_mean_amp] = sort(mean_amp,'ascend'); 

mean_alpha_amp = squeeze(mean(mean(mampdiff(:,postchans,4:12),2),3)); 
[sv_mean_alpha_amp,si_mean_alpha_amp] = sort(mean_alpha_amp,'ascend'); 

snr_gamma = squeeze(mean(mean(mean(allmersp(:,:,20:40,40:80),2),3),4)); 
snr_alpha = squeeze(mean(mean(mean(allmersp(:,:,5:12,40:80),2),3),4)); 

right_chans = [58,26,59,27,63,64,31];
left_chans = [23,56,24,57,60,61,29]; 

mean_sbothdists = squeeze(mean(sbothdists(:,postchans),2)); 

colors = colormap(parula); colors = imresize(colors,[24,3]); colors(colors>1) = 1; colors(colors<0) = 0; 
close all; 
goodsubs = [3,4,5,6,7,8,10,11,12,15,18,19,20,21,22,23];
for i=1:24
    subplot(4,6,i) ; 
    imagesc(squeeze(mean(mean(allmersp(i,[1,3,5],:,:),1),2)),[-1,1]); axis xy;
    erspi = squeeze(mean(mean(mean(allmersp(i,[1,3,5],:,times>0 & times<2),1),2),4)); 
    maxpks(i) = find(erspi(20:end) == max(erspi(20:end))); 
    max_apks(i) = find(erspi(4:12) == min(erspi(4:12))); 
    title(['gp = ',num2str((maxpks(i)+19)*2),', ap = ',num2str((max_apks(i)+3)*2)]);
end

[c,p] = corr(mean(sbothdists(:,postchans),2),mean(inots(:,1:end/2),2)); 

plot(squeeze(mean(mean(mampdiff(:,postchans,4:7),2),3)),squeeze(mean(mean(mampdiff(:,postchans,8:12),2),3)),'kd','LineWidth',2); lsline; 
[c,p] = corr(squeeze(mean(mean(mampdiff(:,postchans,4:7),2),3)),squeeze(mean(mean(mampdiff(:,postchans,8:12),2),3)));



% for peak frequency:
mdist = squeeze(mean(sbothdists(:,postchans),2)); 

[c,p] = corr(mdist(peaksubs),allpeaks(peaksubs)'); 
[c,p] = corr(binots(peaksubs)',allpeaks(peaksubs)'); 


mmersp = squeeze(mean(mean(allmersp(:,[1,5],:,50:80),2),4)); 
for i=1:24 ; minalphas(i) = find(mmersp(i,3:13)==min(mmersp(i,3:13))) ; end
[c,p] = corr(mdist,minalphas'); 
[c,p] = corr(binots(:),minalphas'); 


% figure 1 (methods + data) 
figure,
% plot ersp
subplot(2,3,1);
imagesc(times,freqs,squeeze(mean(mean(allmersp(:,1,:,:),1),2)),[-1,1]); axis xy; 
vline([0,2],'k'); xlabel('time(s)'); ylabel('frequency(hz)'); 

mtopo_alpha = squeeze(mean(mean(mean(stimampdiff(:,1,:,4:12),1),2),4)); 
mtopo_gamma = squeeze(mean(mean(mean(stimampdiff(:,1,:,20:end),1),2),4));

% plot
subplot(2,3,2); 
topoplot(squeeze(mean(mean(mean(stimampdiff(:,1,:,20:end),1),2),4)),eeg.chanlocs,'emarker2',{postchans,'.','r'}); colormap parula;

subplot(2,3,3); 
topoplot(squeeze(mean(mean(mean(stimampdiff(:,1,:,4:12),1),2),4)),eeg.chanlocs,'emarker2',{postchans,'.','r'}); colormap parula;

%for i=1:24
%    subplot(2,12,i); 
%    imagesc(times,freqs,squeeze(mean(allmersp(i,[1],:,:),2)),[-2,2]);  axis xy; colormap parula
%    title(['subject ',num2str(i)]);
%end


subplot(5,3,15); 
bar(vec2(1:50),'b'); hold on ; plot(mat2gray(mean(resdenoise(tcorrinds,1:50),1)),'r','LineWidth',2)
xlabel('time in TR (TR=2s)'); ylabel('BOLD amp. (A.U)'); 

subplot(5,3,14); imagesc([-1,1]) ; colorbar;
subplot(5,3,13); imagesc([min(mtopo_alpha),abs(min(mtopo_alpha))]) ; colorbar; 
subplot(5,3,10); imagesc([-max(mtopo_gamma),max(mtopo_gamma)]) ; colorbar; 

figure,
subplot(2,3,1); 
shadedErrorBar(1:2:120,squeeze(mean(mean(allmersp(:,1,:,times>0 & times<2),1),4)),squeeze(std(mean(allmersp(:,1,:,times>0 & times<2),4),0,1))./sqrt(24),{'Color',[.6,0,0]}); hold on; 
shadedErrorBar(1:2:120,squeeze(mean(mean(allmersp(:,3,:,times>0 & times<2),1),4)),squeeze(std(mean(allmersp(:,3,:,times>0 & times<2),4),0,1))./sqrt(24),{'Color',[0,0,.6]});  
shadedErrorBar(1:2:120,squeeze(mean(mean(allmersp(:,2,:,times>0 & times<2),1),4)),squeeze(std(mean(allmersp(:,2,:,times>0 & times<2),4),0,1))./sqrt(24),{'Color',[0,.6,0]}); hline(0,'k'); 
xlabel('frequency(hz)'); ylabel('db'); xlim([0,120]);ylim([-2.5,1]);

subplot(2,3,2); 
shadedErrorBar(1:2:120,squeeze(mean(mean(allmersp(:,5,:,times>0 & times<2),1),4)),squeeze(std(mean(allmersp(:,1,:,times>0 & times<2),4),0,1))./sqrt(24),{'Color',[.6,0,.6]}); hold on; 
shadedErrorBar(1:2:120,squeeze(mean(mean(allmersp(:,6,:,times>0 & times<2),1),4)),squeeze(std(mean(allmersp(:,3,:,times>0 & times<2),4),0,1))./sqrt(24),{'Color',[0,.6,.6]});  
shadedErrorBar(1:2:120,squeeze(mean(mean(allmersp(:,4,:,times>0 & times<2),1),4)),squeeze(std(mean(allmersp(:,2,:,times>0 & times<2),4),0,1))./sqrt(24),{'Color',[.6,.3,0]}); hline(0,'k'); 
xlabel('frequency(hz)'); ylabel('db'); xlim([0,120]); ylim([-2.5,1]);

subplot(2,3,4); 
plot(1,'Color',[.6,0,0]); hold on; plot(1,'Color',[0,0,.6]); plot(1,'Color',[0,.6,0]); legend({'100% random','60% random','plaid'});
subplot(2,3,5); 
plot(1,'Color',[.6,0,0.6]); hold on; plot(1,'Color',[0,0.6,.6]); plot(1,'Color',[0.6,.3,0]); legend({'10% random','60% random','plaid'});



% figure 2 (inter-frequency and inter-stimulus correlations)
figure,
subplot(3,7,1); 
[mampsv,mampsi] = sort(mean(mean(mean(cleanamps(:,postchans,20:50,:),2),3),4)); 
gamma_amps = mean(mean(mean(cleanamps(:,postchans,20:50,:),2),3),4); 
hold on
for i=1:length(mampsv)
    h=bar(i,mampsv(i));
    ind = find(sort(mampsv)==mampsv(i));     
    set(h,'FaceColor',[colors(ind,1),colors(ind,2),colors(ind,3)]);    
end
errorbar(mampsv,squeeze(std(mean(mean(cleanamps(:,postchans,20:50,:),2),3),0,4))/sqrt(720),'k.'); 
xlabel('subject'); ylabel('gamma amplitude (micro-volts)'); 
title('sorted gamma amp.'); 

subplot(3,7,2);
plot(squeeze(mean(mean(mampdiff(:,postchans,4:7),2),3)),squeeze(mean(mean(mampdiff(:,postchans,8:12),2),3)),'kd','LineWidth',2);  lsline
xlabel('alpha amplitude (micro-volts)'); ylabel('beta amplitude (micro-volts)'); 
[c,p] = corr(squeeze(mean(mean(mampdiff(:,postchans,4:7),2),3)),squeeze(mean(mean(mampdiff(:,postchans,8:12),2),3))); 
[xl,yl] = getlimits(squeeze(mean(mean(mampdiff(:,postchans,4:7),2),3)),squeeze(mean(mean(mampdiff(:,postchans,8:12),2),3))); 
xlim(xl); ylim(yl); 
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]); 

%{
alpha_amps = (mean(mean(mean(cleanamps(:,postchans,4:12,:),2),3),4)); 
hold on
for i=1:length(alpha_amps)
    h=bar(i,alpha_amps(mampsi(i)));
    ind = find(sort(alpha_amps)==alpha_amps(mampsi(i)));     
    set(h,'FaceColor',[colors(ind,1),colors(ind,2),colors(ind,3)]);    
end
errorbar(alpha_amps(mampsi),squeeze(std(mean(mean(cleanamps(mampsi,postchans,4:12,:),2),3),0,4))/sqrt(720),'k.'); 
xlabel('subject'); ylabel('alpha/beta amplitude (micro-volts)'); 
title('alpha/beta amp. sorted by gamma amp.'); 
%}

subplot(3,7,3); 
plot(gamma_amps,alpha_amps,'kd','LineWidth',2);  lsline
xlabel('gamma amplitude (micro-volts)'); ylabel('alpha/beta amplitude (micro-volts)'); 
[c,p] = corr(gamma_amps,alpha_amps); 
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]); 

subplot(3,7,4); 
titles = {'100% contrast','5% contrast','33% contrast','plaid','10% random','60% random'};
pointforms = {'rd','bd','md','kd','gd'};lineforms = {'-r','-b','-m','-k','-g'};
clear rhos ps; 
for i=2:6
x1 = squeeze(mean(mean(smampdiff(:,1,postchans,20:end),3),4)); y1 = squeeze(mean(mean(smampdiff(:,i,postchans,20:end),3),4)); 
p = polyfit(x1,y1,1) ; r = p(1) .* x1 + p(2) ; 
plot(x1,y1,pointforms{i-1},'LineWidth',2); hold on; 
plot(x1,r,lineforms{i-1},'LineWidth',2) ; 
[rhos(i-1),ps(i-1)] = corr(x1,y1); 
end
xlabel('gamma amplitude'); ylabel('gamma amplitude'); xlim([-0.002,0.02]);
title(['median rho=',num2str(round(median(rhos*100))/100),', median ',format_p(median(ps))]); 

subplot(3,7,5); 
pointforms = {'rd','bd','md','kd','gd'};lineforms = {'-r','-b','-m','-k','-g'};
clear rhos ps; 
for i=2:6
x1 = squeeze(mean(mean(smampdiff(:,1,postchans,4:12),3),4)); y1 = squeeze(mean(mean(smampdiff(:,i,postchans,4:12),3),4)); 
p = polyfit(x1,y1,1) ; r = p(1) .* x1 + p(2) ; 
plot(x1,y1,pointforms{i-1},'LineWidth',2); hold on; 
plot(x1,r,lineforms{i-1},'LineWidth',2) ; 
[rhos(i-1),ps(i-1)] = corr(x1,y1); 
end
xlabel('alpha/beta amplitude'); ylabel('alpha/beta amplitude'); 
title(['median rho=',num2str(round(median(rhos*100))/100),', median ',format_p(median(ps))]); 

subplot(3,7,6) ; plot(1,'r'); hold on ; plot(1,'b'); plot(1,'m'); plot(1,'k'); plot(1,'g'); 
legend({'vs 5% contrast','vs 33% contrast','vs plaid','vs 10% random','vs 60% random'});

[rhos,ps] = corr(squeeze(mean(mean(smampdiff(:,:,postchans,:),2),3))); 
subplot(3,7,7); imagesc(1:2:100,1:2:100,rhos,[-1,1]) ; xlabel('frequency(hz)'); ylabel('frequency(hz)'); 

stimindices = [1,3,2,5,6,4]; 
for i=1:6 ; subplot(5,14,42+i) ; imagesc(times,freqs,squeeze(allmersp(21,stimindices(i),:,:)),[-2,2]); axis xy; if i==1; xlabel('time(s)'); ylabel('frequency(hz)');end; title(titles{stimindices(i)}); end 
for i=1:6 ; subplot(5,14,56+i) ; imagesc(times,freqs,squeeze(allmersp(16,stimindices(i),:,:)),[-2,2]); axis xy; end;  subplot(4,12,46) ; imagesc([-2,2]) ; colorbar; 
subplot(4,12,47) ; imagesc([-1,1]); colorbar;

% badsubs
subplot(2,2,1);
plot(squeeze(mean(sbothdists(goodsubs,postchans),2)),squeeze(mean(mean(mampdiff(goodsubs,postchans,20:50),2),3)),'kd','LineWidth',2); lsline
[c,p] = corr(squeeze(mean(sbothdists(goodsubs,postchans),2)),squeeze(mean(mean(mampdiff(goodsubs,postchans,20:50),2),3))); 
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]); xlabel('distance (mm)'); ylabel('gamma amplitude (micro-volts)'); 
[xl,yl] = getlimits(squeeze(mean(sbothdists(goodsubs,postchans),2)),squeeze(mean(mean(mampdiff(goodsubs,postchans,20:50),2),3))); 
xlim(xl); ylim(yl); 



%{
% figure 3 (distance vs eeg) 
figure,

subplot(2,2,1);
plot(squeeze(mean(sbothdists(:,postchans),2)),squeeze(mean(mean(mampdiff(:,postchans,20:40),2),3)),'kd','LineWidth',2); lsline
[c,p] = corr(squeeze(mean(sbothdists(:,postchans),2)),squeeze(mean(mean(mampdiff(:,postchans,20:40),2),3))); 
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]); xlabel('distance (mm)'); ylabel('gamma amplitude (micro-volts)'); 
[xl,yl] = getlimits(squeeze(mean(sbothdists(:,postchans),2)),squeeze(mean(mean(mampdiff(:,postchans,20:40),2),3))); 
xlim(xl); ylim(yl); 

subplot(2,2,2); 
plot(squeeze(mean(sbothdists(:,postchans),2)),squeeze(mean(mean(mampdiff(:,postchans,5:12),2),3)),'kd','LineWidth',2); lsline
[c,p] = corr(squeeze(mean(sbothdists(:,postchans),2)),squeeze(mean(mean(mampdiff(:,postchans,5:12),2),3))); 
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]); xlabel('distance (mm)'); ylabel('alpha/beta amplitude (micro-volts)'); 
[xl,yl] = getlimits(squeeze(mean(sbothdists(:,postchans),2)),squeeze(mean(mean(mampdiff(:,postchans,5:12),2),3))); 
xlim(xl); ylim(yl); 

subplot(4,4,9);
[dist_sv,dist_si] = sort(mean_sbothdists,'descend'); 
hold on
for i=1:length(sv_mean_amp)
    h=bar(i,mean_amp(dist_si(i)));
    ind = find(sort(sv_mean_amp)==mean_amp(dist_si(i)));     
    set(h,'FaceColor',[colors(ind,1),colors(ind,2),colors(ind,3)]);    
end
xlim([0,25]); xlabel('subject'); ylabel('gamma amplitude'); 
title('gamma amp. sorted by distance'); box on;
[xl,yl] = getlimits(ones(1,10),mean_amp); ylim(yl); 

subplot(4,4,10);
group_dist_gamma = [squeeze((mean_amp(dist_si(1:11)))),squeeze((mean_amp(dist_si(end-10:end))))];
bar(squeeze(mean(group_dist_gamma,1))); hold on;
h=bar(1,mean(group_dist_gamma(:,1))); set(h,'FaceColor',[colors(1,1),colors(1,2),colors(1,3)]); 
h=bar(2,mean(group_dist_gamma(:,2))); set(h,'FaceColor',[colors(end,1),colors(end,2),colors(end,3)]); 
errorbar(squeeze(mean(group_dist_gamma,1)),squeeze(std(group_dist_gamma,0,1))./sqrt(24),'k.'); 
group_distlabs = [min(mean_sbothdists(dist_si(1:11))),max(mean_sbothdists(dist_si(1:11))),min(mean_sbothdists(dist_si(end-10:end))),max(mean_sbothdists(dist_si(end-10:end)))]; 
xlim([0.5,2.5]); xlabel('distance group'); ylabel('gamma amplitude'); set(gca,'XTickLabel',{'high(56-61mm)','low(51-55mm)'}); 
[h,p,ci,stats] = ttest(group_dist_gamma(:,1),group_dist_gamma(:,2)); 
title(['t=',num2str(round(stats.tstat*100)/100),', ',format_p(p)]);

subplot(4,4,11); 
[dist_sv,dist_si] = sort(mean_sbothdists,'descend'); 
hold on
for i=1:length(sv_mean_alpha_amp)
    h=bar(i,mean_alpha_amp(dist_si(i)));
    ind = find(sort(sv_mean_alpha_amp)==mean_alpha_amp(dist_si(i)));     
    set(h,'FaceColor',[colors(ind,1),colors(ind,2),colors(ind,3)]);    
end
xlim([0,25]); xlabel('subject'); ylabel('alpha/beta amplitude'); 
title('alpha/beta amp. sorted by distance'); box on;

subplot(4,4,12); 
group_dist_alpha = [squeeze((mean_alpha_amp(dist_si(1:11)))),squeeze((mean_alpha_amp(dist_si(end-10:end))))];
bar(squeeze(mean(group_dist_alpha,1))); hold on;
h=bar(1,mean(group_dist_alpha(:,1))); set(h,'FaceColor',[colors(end,1),colors(end,2),colors(end,3)]); 
h=bar(2,mean(group_dist_alpha(:,2))); set(h,'FaceColor',[colors(1,1),colors(1,2),colors(1,3)]); 
errorbar(squeeze(mean(group_dist_alpha,1)),squeeze(std(group_dist_alpha,0,1))./sqrt(24),'k.'); 
xlim([0.5,2.5]); xlabel('distance group'); ylabel('alpha/beta amplitude'); set(gca,'XTickLabel',{'high(56-61mm)','low(51-55mm)'});  
[h,p,ci,stats] = ttest(group_dist_alpha(:,1),group_dist_alpha(:,2)); 
title(['t=',num2str(round(stats.tstat*100)/100),', ',format_p(p)]);

subplot(4,4,13);
topoplot(gamma_ecorrs,eeg.chanlocs,'maplimits',[-.7,.7],'style','both','emarker2',{find(gamma_eps(goodchans)<0.05),'.','w'},'plotchans',goodchans); colormap parula;
title('40-100Hz'); 

subplot(4,4,14); 
hemi_elecdists = [squeeze(mean(srightdists(:,right_chans),2)),squeeze(mean(sleftdists(:,left_chans),2))]; 
bar(squeeze(mean(hemi_elecdists,1)),'c'); hold on; 
errorbar(squeeze(mean(hemi_elecdists,1)),squeeze(std(hemi_elecdists,0,1))./sqrt(24),'k.'); 
ylim([40,55]); xlim([0.5,2.5]);
[h,p,ci,stats] = ttest(hemi_elecdists(:,1),hemi_elecdists(:,2)); 
title(['t=',num2str(round(stats.tstat*100)/100),', ',format_p(p)]); set(gca,'XTickLabel',{'right','left'}); ylabel('distance(mm)'); xlabel('hemisphere');

subplot(4,4,15); 
errorbarxy(squeeze(mean(sbothdists(:,goodchans),1)),squeeze(mean(mean((mampdiff(:,goodchans,20:50)),1),3)),...
squeeze(std(sbothdists(:,goodchans),0,1))./sqrt(24),squeeze(std(mean((mampdiff(:,goodchans,20:50)),3),0,1))./sqrt(24),{'kd','k','k'});
ylim([0,0.007]); xlim([20,180]); xlabel('distance(mm)'); ylabel('gamma amp. (microvolts)');

subplot(4,4,16); 
errorbarxy(squeeze(mean(sbothdists(:,goodchans),1)),squeeze(mean(mean((mampdiff(:,goodchans,4:12)),1),3)),...
squeeze(std(sbothdists(:,goodchans),0,1))./sqrt(24),squeeze(std(mean((mampdiff(:,goodchans,4:12)),3),0,1))./sqrt(24),{'kd','k','k'});
ylim([-0.32,0]); xlim([20,180]); xlabel('distance(mm)'); ylabel('alpha/beta amp. (microvolts)');  

%}


% figure 4 (i0 vs gamma and alpha/beta)
[icorrs,ips] = corr(inots,squeeze(mean(mampdiff(:,postchans,:),2))); 

subplot(2,2,1); 
plot(squeeze(mean(inots(:,1:end/2),2)),squeeze(mean(mean(mampdiff(:,postchans,20:50),2),3)),'kd','LineWidth',2); lsline
[c,p] = corr(squeeze(mean(inots(:,1:end/2),2)),squeeze(mean(mean(mampdiff(:,postchans,20:50),2),3))); 
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]); 
xlabel('I0'); ylabel('gamma amplitude (micro-volts)'); 
[xl,yl] = getlimits(squeeze(mean(inots(:,1:end/2),2)),squeeze(mean(mean(mampdiff(:,postchans,20:50),2),3))); 
xlim(xl); ylim(yl); 

subplot(2,2,2); 
plot(squeeze(mean(inots(:,1:end/2),2)),squeeze(mean(mean(mampdiff(:,postchans,4:12),2),3)),'kd','LineWidth',2); lsline
[c,p] = corr(squeeze(mean(inots(:,1:end/2),2)),squeeze(mean(mean(mampdiff(:,postchans,4:12),2),3))); 
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]); 
xlabel('I0'); ylabel('alpha/beta amplitude (micro-volts)'); 
[xl,yl] = getlimits(squeeze(mean(inots(:,1:end/2),2)),squeeze(mean(mean(mampdiff(:,postchans,4:12),2),3))); 
xlim(xl); ylim(yl); 

% bar distance sorted gamma
subplot(4,4,9) 
[inot_sv,inot_si] = sort(squeeze(mean(inots(:,1:end/2),2)),'descend'); 
hold on
for i=1:length(sv_mean_amp)
    h=bar(i,mean_amp(inot_si(i)));
    ind = find(sort(sv_mean_amp)==mean_amp(inot_si(i)));     
    set(h,'FaceColor',[colors(ind,1),colors(ind,2),colors(ind,3)]);    
end
xlim([0,25]); xlabel('subject'); ylabel('gamma amplitude'); 
title('gamma amp. sorted by I0'); box on;
[xl,yl] = getlimits(ones(1,10),mean_amp); ylim(yl); 
% bar group difference distance sorted gamma
subplot(4,4,10) 
group_i0_gamma = [squeeze((mean_amp(inot_si(1:11)))),squeeze((mean_amp(inot_si(end-10:end))))];
bar(squeeze(mean(group_i0_gamma,1))); hold on;
h=bar(1,mean(group_i0_gamma(:,1))); set(h,'FaceColor',[colors(1,1),colors(1,2),colors(1,3)]); 
h=bar(2,mean(group_i0_gamma(:,2))); set(h,'FaceColor',[colors(end,1),colors(end,2),colors(end,3)]); 
errorbar(squeeze(mean(group_i0_gamma,1)),squeeze(std(group_i0_gamma,0,1))./sqrt(24),'k.'); 
group_distlabs = [min(mean_sbothdists(dist_si(1:11))),max(mean_sbothdists(dist_si(1:11))),min(mean_sbothdists(dist_si(end-10:end))),max(mean_sbothdists(dist_si(end-10:end)))]; 
xlim([0.5,2.5]); xlabel('I0 group'); ylabel('gamma amplitude'); set(gca,'XTickLabel',{'high(0.95-1)','low(0.9-0.94)'}); 
[h,p,ci,stats] = ttest(group_i0_gamma(:,1),group_i0_gamma(:,2)); 
title(['t=',num2str(round(stats.tstat*100)/100),', ',format_p(p)]);


% bar i0 sorted alpha/beta
subplot(4,4,11) 
[inot_sv,inot_si] = sort(squeeze(mean(inots(:,1:end/2),2)),'descend'); 
hold on
for i=1:length(sv_mean_alpha_amp)
    h=bar(i,mean_alpha_amp(inot_si(i)));
    ind = find(sort(sv_mean_alpha_amp)==mean_alpha_amp(inot_si(i)));     
    set(h,'FaceColor',[colors(ind,1),colors(ind,2),colors(ind,3)]);    
end
xlim([0,25]); xlabel('subject'); ylabel('alpha amplitude'); 
title('alpha amp. sorted by I0'); box on;
[xl,yl] = getlimits(ones(1,10),mean_alpha_amp); ylim(yl); 
% bar group difference i0 sorted alpha/beta
subplot(4,4,12) 
group_i0_alpha = [squeeze((mean_alpha_amp(inot_si(1:11)))),squeeze((mean_alpha_amp(inot_si(end-10:end))))];
bar(squeeze(mean(group_i0_alpha,1))); hold on;
h=bar(1,mean(group_i0_alpha(:,1))); set(h,'FaceColor',[colors(end,1),colors(end,2),colors(end,3)]); 
h=bar(2,mean(group_i0_alpha(:,2))); set(h,'FaceColor',[colors(1,1),colors(1,2),colors(1,3)]); 
errorbar(squeeze(mean(group_i0_alpha,1)),squeeze(std(group_i0_alpha,0,1))./sqrt(24),'k.'); 
group_distlabs = [min(mean_sbothdists(dist_si(1:11))),max(mean_sbothdists(dist_si(1:11))),min(mean_sbothdists(dist_si(end-10:end))),max(mean_sbothdists(dist_si(end-10:end)))]; 
xlim([0.5,2.5]); xlabel('I0 group'); ylabel('alpha amplitude'); set(gca,'XTickLabel',{'high(0.95-1)','low(0.9-0.94)'}); 
[h,p,ci,stats] = ttest(group_i0_alpha(:,1),group_i0_alpha(:,2)); 
title(['t=',num2str(round(stats.tstat*100)/100),', ',format_p(p)]);


[alpha_inotcorrs,alpha_inotps] = corr(mean(inots(:,1:end/2),2),mean(mampdiff(:,:,4:12),3)); 

%{
subplot(4,4,11);
topoplot(alpha_inotcorrs,eeg.chanlocs,'maplimits',[-.7,.7],'style','both','emarker2',{find(alpha_inotps(goodchans)<0.05),'.','w'},'plotchans',goodchans,'maplimits',[-.5,.5]); colormap parula;
title('8-25'); 

subplot(3,5,9); 
hemi_elecdists = [squeeze(mean(srightdists(:,right_chans),2)),squeeze(mean(sleftdists(:,left_chans),2))]; 
bar(squeeze(mean(hemi_elecdists,1)),'c'); hold on; 
errorbar(squeeze(mean(hemi_elecdists,1)),squeeze(std(hemi_elecdists,0,1))./sqrt(24),'k.'); 
ylim([40,55]); xlim([0.5,2.5]);
[h,p,ci,stats] = ttest(hemi_elecdists(:,1),hemi_elecdists(:,2)); 
title(['t=',num2str(round(stats.tstat*100)/100),', ',format_p(p)]); set(gca,'XTickLabel',{'right','left'}); ylabel('distance(mm)'); xlabel('hemisphere');

subplot(3,5,10);
shadedErrorBar([],squeeze(mean(inots,1)),squeeze(std(inots,0,1))/sqrt(24)); 
xlim([1,50]);
%}


















%{
figure,
% corr alpha/beta amp vs dist
subplot(3,11,1);
plot(squeeze(mean(sbothdists(:,postchans),2)),squeeze(mean(mean(mampdiff(:,postchans,5:12),2),3)),'kd','LineWidth',1); lsline
[c,p] = corr(squeeze(mean(sbothdists(:,postchans),2)),squeeze(mean(mean(mampdiff(:,postchans,5:12),2),3))); 
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]); xlabel('distance (mm)'); ylabel('amplitude (micro-volts)'); 
[xl,yl] = getlimits(squeeze(mean(sbothdists(:,postchans),2)),squeeze(mean(mean(mampdiff(:,postchans,5:12),2),3))); 
xlim(xl); ylim(yl); 
% corr gamma amp vs dist
subplot(3,11,2);
plot(squeeze(mean(sbothdists(:,postchans),2)),squeeze(mean(mean(mampdiff(:,postchans,20:40),2),3)),'kd','LineWidth',1); lsline
[c,p] = corr(squeeze(mean(sbothdists(:,postchans),2)),squeeze(mean(mean(mampdiff(:,postchans,20:40),2),3))); 
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]); xlabel('distance (mm)'); ylabel('amplitude (micro-volts)'); 
[xl,yl] = getlimits(squeeze(mean(sbothdists(:,postchans),2)),squeeze(mean(mean(mampdiff(:,postchans,20:40),2),3))); 
xlim(xl); ylim(yl); 
% corr alpha/beta SNR vs distance
subplot(3,11,3);
plot(squeeze(mean(sbothdists(:,postchans),2)),snr_alpha,'kd','LineWidth',1); lsline
[c,p] = corr(squeeze(mean(sbothdists(:,postchans),2)),snr_alpha); 
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]); xlabel('distance (mm)'); ylabel('SNR (db)'); 
[xl,yl] = getlimits(squeeze(mean(sbothdists(:,postchans),2)),snr_alpha); 
xlim(xl); ylim(yl); 
% corr gamma SNR vs distance
subplot(3,11,4);
plot(squeeze(mean(sbothdists(:,postchans),2)),snr_gamma,'kd','LineWidth',1); lsline
[c,p] = corr(squeeze(mean(sbothdists(:,postchans),2)),snr_gamma); 
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]); xlabel('distance (mm)'); ylabel('SNR (db)'); 
[xl,yl] = getlimits(squeeze(mean(sbothdists(:,postchans),2)),snr_gamma); 
xlim(xl); ylim(yl); 
% corr alpha amp vs gamma amp
subplot(3,11,5);
plot(squeeze(mean(mean(mampdiff(:,postchans,20:40),2),3)),squeeze(mean(mean(mampdiff(:,postchans,5:12),2),3)),'kd','LineWidth',1); lsline
[c,p] = corr(squeeze(mean(mean(mampdiff(:,postchans,20:40),2),3)),squeeze(mean(mean(mampdiff(:,postchans,5:12),2),3))); 
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]); xlabel('amplitude (micro-volts)'); ylabel('amplitude (micro-volts)'); 
[xl,yl] = getlimits(squeeze(mean(mean(mampdiff(:,postchans,20:40),2),3)),squeeze(mean(mean(mampdiff(:,postchans,5:12),2),3))); 
xlim(xl); ylim(yl); 
% corr alpha SNR vs gamma SNR
subplot(3,11,6);
plot(snr_alpha,snr_gamma,'kd','LineWidth',1); lsline
[c,p] = corr(snr_alpha,snr_gamma); 
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]); xlabel('amplitude (micro-volts)'); ylabel('amplitude (micro-volts)'); 
[xl,yl] = getlimits(snr_alpha,snr_gamma); 
xlim(xl); ylim(yl); 
% bar amplitude gamma vs amplitude alpha 
subplot(3,11,7);
m_gamma = squeeze(mean(mean(mampdiff(:,postchans,20:40),2),3)); 
m_alpha = squeeze(mean(mean(mampdiff(:,postchans,5:12),2),3)); 
barwitherr(std([m_gamma,m_alpha],0,1)/sqrt(24),mean([m_gamma,m_alpha],1),'c'); 
amp_ratio = mean(m_alpha)./mean(m_gamma); 
title(['amplitude ratio = ',num2str(abs(round(amp_ratio*100)/100))]);  set(gca,'XTickLabel',{'40-90','8-25'}); xlabel('frequency band (hz)'); ylabel('amplitude (micro-volts)'); 
xlim([0.5,2.5]); ylim([-0.3,0.05]);
% bar SNR gamma vs alpha
subplot(3,11,8);
barwitherr(std([snr_gamma,snr_alpha],0,1)/sqrt(24),mean([snr_gamma,snr_alpha],1),'c'); 
snr_ratio = mean(snr_alpha)./mean(snr_gamma); 
title(['SNR ratio = ',num2str(abs(round(snr_ratio*100)/100))]); set(gca,'XTickLabel',{'40-90','8-25'}); xlabel('frequency band (hz)'); ylabel('SNR (db)'); 
xlim([0.5,2.5]); ylim([-1.5,0.5]);
% plot snr vs alpha amplitude
subplot(3,11,9);
plot(snr_alpha,m_alpha,'kd','LineWidth',1); lsline; 
[c,p] = corr(m_alpha,snr_alpha); 
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]); xlabel('SNR (db)'); ylabel('amplitude (micro-volt)'); 
[xl,yl] = getlimits(snr_alpha,m_alpha); 
xlim(xl); ylim(yl); 
% plot SNR vs gamma amplitude
subplot(3,11,10);
plot(snr_gamma,m_gamma,'kd','LineWidth',1); lsline; 
[c,p] = corr(m_gamma,snr_gamma); 
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]); xlabel('SNR (db)'); ylabel('amplitude (micro-volt)'); 
[xl,yl] = getlimits(snr_gamma,m_gamma); 
xlim(xl); ylim(yl); 
% plot SNR from 1-120Hz
subplot(3,11,11);
merspfreqs = squeeze(mean(mean(allmersp(:,:,:,40:80),2),4)); 
shadedErrorBar(1:2:120,mean(abs(merspfreqs)),std(abs(merspfreqs))/sqrt(24)) 
ylim([-0.5,2.5]); hline(0,'k'); title('SNR 1-120Hz'); xlim([0,121]);  xlabel('frequency (hz)'); ylabel('SNR (db)'); 
% plot amplitude from 1-120Hz
subplot(3,11,12); 
shadedErrorBar(1:2:100,squeeze(mean(mean(abs(mampdiff(:,postchans,:)),1),2)),squeeze(std(mean(abs(mampdiff(:,postchans,:)),2),0,1))./sqrt(24)) 
title('amplitude 1-100Hz'); xlabel('frequency(hz)'); ylabel('amplitude (micro-volts)'); ylim([0,1]); 
% plot alpha amp as a function of distance (all electrodes)
subplot(3,11,13); 
errorbarxy(squeeze(mean(sbothdists(:,goodchans),1)),squeeze(mean(mean((mampdiff(:,goodchans,4:12)),1),3)),...
squeeze(std(sbothdists(:,goodchans),0,1))./sqrt(24),squeeze(std(mean((mampdiff(:,goodchans,4:12)),3),0,1))./sqrt(24),{'kd','k','k'});
ylim([-0.32,0]); xlim([20,180]); title('8-25Hz'); xlabel('distance(mm)'); ylabel('alpha amp. (microvolts)');  
% plot gamma amp as a function of distance (all electrodes)
subplot(3,11,14); 
errorbarxy(squeeze(mean(sbothdists(:,goodchans),1)),squeeze(mean(mean((mampdiff(:,goodchans,20:50)),1),3)),...
squeeze(std(sbothdists(:,goodchans),0,1))./sqrt(24),squeeze(std(mean((mampdiff(:,goodchans,20:50)),3),0,1))./sqrt(24),{'kd','k','k'});
ylim([0,0.007]); xlim([20,180]); title('40-100Hz');xlabel('distance(mm)'); ylabel('gamma amp. (microvolts)');
% topoplot alpha rho 
subplot(3,11,15); 
topoplot(alpha_ecorrs,eeg.chanlocs,'maplimits',[-.7,.7],'style','both','plotchans',goodchans); colormap parula;
title('8-25Hz');
% topoplot gamma rho 
subplot(3,11,16); 
topoplot(gamma_ecorrs,eeg.chanlocs,'maplimits',[-.7,.7],'style','both','emarker2',{find(gamma_eps(goodchans)<0.05),'.','w'},'plotchans',goodchans); colormap parula;
title('40-100Hz'); 
% bar mean distance left and right
subplot(3,11,17); 
hemi_elecdists = [squeeze(mean(srightdists(:,right_chans),2)),squeeze(mean(sleftdists(:,left_chans),2))]; 
bar(squeeze(mean(hemi_elecdists,1)),'c'); hold on; 
errorbar(squeeze(mean(hemi_elecdists,1)),squeeze(std(hemi_elecdists,0,1))./sqrt(24),'k.'); 
ylim([40,55]); xlim([0.5,2.5]);
[h,p,ci,stats] = ttest(hemi_elecdists(:,1),hemi_elecdists(:,2)); 
title(['t=',num2str(round(stats.tstat*100)/100),', ',format_p(p)]); set(gca,'XTickLabel',{'right','left'}); ylabel('distance(mm)'); xlabel('hemisphere');
% bar mean alpha power left and right
subplot(3,11,18)
hemi_amp = [squeeze(mean(mean(mampdiff(:,right_chans,4:12),2),3)),squeeze(mean(mean(mampdiff(:,left_chans,4:12),2),3))];
bar(squeeze(mean(hemi_amp,1)),'c'); hold on; 
errorbar(squeeze(mean(hemi_amp,1)),squeeze(std(hemi_amp,0,1))./sqrt(24),'k.'); 
[h,p,ci,stats] = ttest(hemi_amp(:,1),hemi_amp(:,2)); 
title(['t=',num2str(round(stats.tstat*100)/100),', ',format_p(p)]); set(gca,'XTickLabel',{'right','left'}); ylabel('distance(mm)'); xlabel('hemisphere');
% bar mean gamma power left and right
subplot(3,11,19) 
hemi_amp = [squeeze(mean(mean(mampdiff(:,right_chans,20:50),2),3)),squeeze(mean(mean(mampdiff(:,left_chans,20:50),2),3))];
bar(squeeze(mean(hemi_amp,1)),'c'); hold on; 
errorbar(squeeze(mean(hemi_amp,1)),squeeze(std(hemi_amp,0,1))./sqrt(24),'k.'); 
[h,p,ci,stats] = ttest(hemi_amp(:,1),hemi_amp(:,2)); 
title(['t=',num2str(round(stats.tstat*100)/100),', ',format_p(p)]); set(gca,'XTickLabel',{'right','left'}); ylabel('distance(mm)'); xlabel('hemisphere');

% bar gamma sorted gamma
subplot(3,11,20) 
hold on
for i=1:length(sv_mean_amp)
    h=bar(i,sv_mean_amp(i));
    ind = find(sort(sv_mean_amp)==sv_mean_amp(i));     
    set(h,'FaceColor',[colors(ind,1),colors(ind,2),colors(ind,3)]);    
end
xlim([0,25]); ylabel('gamma amplitude'); xlabel('subject'); 
title('sorted gamma amp.'); box on;
[xl,yl] = getlimits(ones(1,10),mean_amp); ylim(yl); 

% bar distance sorted gamma
subplot(3,11,21) 
[dist_sv,dist_si] = sort(mean_sbothdists,'descend'); 
hold on
for i=1:length(sv_mean_amp)
    h=bar(i,mean_amp(dist_si(i)));
    ind = find(sort(sv_mean_amp)==mean_amp(dist_si(i)));     
    set(h,'FaceColor',[colors(ind,1),colors(ind,2),colors(ind,3)]);    
end
xlim([0,25]); xlabel('subject'); ylabel('gamma amplitude'); 
title('gamma amp. sorted by distance'); box on;
[xl,yl] = getlimits(ones(1,10),mean_amp); ylim(yl); 
% bar group difference distance sorted gamma
subplot(3,11,22) 
group_dist_gamma = [squeeze((mean_amp(dist_si(1:11)))),squeeze((mean_amp(dist_si(end-10:end))))];
bar(squeeze(mean(group_dist_gamma,1))); hold on;
h=bar(1,mean(group_dist_gamma(:,1))); set(h,'FaceColor',[colors(1,1),colors(1,2),colors(1,3)]); 
h=bar(2,mean(group_dist_gamma(:,2))); set(h,'FaceColor',[colors(end,1),colors(end,2),colors(end,3)]); 
errorbar(squeeze(mean(group_dist_gamma,1)),squeeze(std(group_dist_gamma,0,1))./sqrt(24),'k.'); 
group_distlabs = [min(mean_sbothdists(dist_si(1:11))),max(mean_sbothdists(dist_si(1:11))),min(mean_sbothdists(dist_si(end-10:end))),max(mean_sbothdists(dist_si(end-10:end)))]; 
xlim([0.5,2.5]); xlabel('distance group'); ylabel('gamma amplitude'); set(gca,'XTickLabel',{'high(56-61mm)','low(51-55mm)'}); 
[h,p,ci,stats] = ttest(group_dist_gamma(:,1),group_dist_gamma(:,2)); 
title(['t=',num2str(round(stats.tstat*100)/100),', ',format_p(p)]);

% bar distance sorted alpha
subplot(3,11,23) 
[dist_sv,dist_si] = sort(mean_sbothdists,'descend'); 
hold on
for i=1:length(sv_mean_alpha_amp)
    h=bar(i,mean_alpha_amp(dist_si(i)));
    ind = find(sort(sv_mean_alpha_amp)==mean_alpha_amp(dist_si(i)));     
    set(h,'FaceColor',[colors(ind,1),colors(ind,2),colors(ind,3)]);    
end
xlim([0,25]); xlabel('subject'); ylabel('alpha amplitude'); 
title('alpha amp. sorted by distance'); box on;

% bar group difference distance sorted alpha
subplot(3,11,24) 
group_dist_alpha = [squeeze((mean_alpha_amp(dist_si(1:11)))),squeeze((mean_alpha_amp(dist_si(end-10:end))))];
bar(squeeze(mean(group_dist_alpha,1))); hold on;
h=bar(1,mean(group_dist_alpha(:,1))); set(h,'FaceColor',[colors(end,1),colors(end,2),colors(end,3)]); 
h=bar(2,mean(group_dist_alpha(:,2))); set(h,'FaceColor',[colors(1,1),colors(1,2),colors(1,3)]); 
errorbar(squeeze(mean(group_dist_alpha,1)),squeeze(std(group_dist_alpha,0,1))./sqrt(24),'k.'); 
xlim([0.5,2.5]); xlabel('distance range(mm)'); ylabel('alpha amplitude'); set(gca,'XTickLabel',{'56-61','51-55'}); 
[h,p,ci,stats] = ttest(group_dist_alpha(:,1),group_dist_alpha(:,2)); 
title(['t=',num2str(round(stats.tstat*100)/100),', ',format_p(p)]);

% bar gamma sorted alpha 
subplot(3,11,25) 
[amp_sv,amp_si] = sort(mean_amp,'descend'); 
hold on
for i=1:length(sv_mean_alpha_amp)
    h=bar(i,mean_alpha_amp(amp_si(i)));
    ind = find(sort(sv_mean_alpha_amp)==mean_alpha_amp(amp_si(i)));     
    set(h,'FaceColor',[colors(ind,1),colors(ind,2),colors(ind,3)]);    
end
xlim([0,25]); xlabel('subject'); ylabel('alpha amplitude'); 
title('alpha amp. sorted by gamma amp.'); box on; 

% bar group difference gamma sorted alpha 
subplot(3,11,26) 
group_dist_alpha = [squeeze((mean_alpha_amp(amp_si(1:11)))),squeeze((mean_alpha_amp(amp_si(end-10:end))))];
bar(squeeze(mean(group_dist_alpha,1))); hold on;
h=bar(1,mean(group_dist_alpha(:,1))); set(h,'FaceColor',[colors(end,1),colors(end,2),colors(end,3)]); 
h=bar(2,mean(group_dist_alpha(:,2))); set(h,'FaceColor',[colors(1,1),colors(1,2),colors(1,3)]); 
errorbar(squeeze(mean(group_dist_alpha,1)),squeeze(std(group_dist_alpha,0,1))./sqrt(24),'k.'); 
group_gammalabs = [min(mean_amp(amp_si(1:11))),max(mean_amp(amp_si(1:11))),min(mean_amp(amp_si(end-10:end))),max(mean_amp(amp_si(end-10:end)))]; 
xlim([0.5,2.5]); xlabel('gamma range(micro volts)'); ylabel('alpha amplitude'); set(gca,'XTickLabel',{'0-0.004','0.005-0.009'}); 
[h,p,ci,stats] = ttest(group_dist_alpha(:,1),group_dist_alpha(:,2)); 
title(['t=',num2str(round(stats.tstat*100)/100),', ',format_p(p)]);
% plot electrodes > 100mm gamma
subplot(3,11,27)
far_elecs = find(mean(sbothdists,1)>100); 
plot(squeeze(mean(sbothdists(:,far_elecs),1)),squeeze(mean(mean(mampdiff(:,far_elecs,4:12),1),3)),'kd'); lsline; 
[c,p] = corr(squeeze(mean(sbothdists(:,far_elecs),1))',squeeze(mean(mean(mampdiff(:,far_elecs,4:12),1),3))');
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]);
% plot electrodes > 100mm alpha/beta
subplot(3,11,28); 
plot(squeeze(mean(sbothdists(:,far_elecs),1)),squeeze(mean(mean(mampdiff(:,far_elecs,20:50),1),3)),'kd'); lsline; 
[c,p] = corr(squeeze(mean(sbothdists(:,far_elecs),1))',squeeze(mean(mean(mampdiff(:,far_elecs,20:50),1),3))');
title(['rho=',num2str(round(c*100)/100),', ',format_p(p)]);


%}


%{
figure,
% plot all ersp 
for i=1:6  
    subplot(3,11,i);    
    imagesc(times,freqs,squeeze(mean(allmersp(:,i,:,:),1)),[-1,1]);  axis xy; 
    if i==1; xlabel('time(s)'); ylabel('frequency(hz)'); end
    subplot(3,11,11+i);
    topoplot(squeeze(mean(mean(stimampdiff(:,i,:,20:40),1),4)),eeg.chanlocs,'maplimits',[-0.005,0.01],'style','map'); colormap parula
    subplot(3,11,22+i); 
    topoplot(squeeze(mean(mean(stimampdiff(:,i,:,4:12),1),4)),eeg.chanlocs,'maplimits',[-0.3,0.01],'style','map'); colormap parula    
end
subplot(3,11,7);    
imagesc(times,freqs,squeeze(mean(mean(allmersp(:,:,:,:),1),2)),[-.5,.5]);  axis xy; 
if i==1; xlabel('time(s)'); ylabel('frequency(hz)'); end
subplot(3,11,11+7);
topoplot(squeeze(mean(mean(mean(stimampdiff(:,:,:,20:50),1),4),2)),eeg.chanlocs,'maplimits',[-0.005,0.006]); colormap parula
subplot(3,11,22+7); 
topoplot(squeeze(mean(mean(mean(stimampdiff(:,:,:,4:12),1),4),2)),eeg.chanlocs,'maplimits',[-0.25,0.01]); colormap parula    
%}



