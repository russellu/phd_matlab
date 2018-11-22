clear all ; close all ; 
cd C:\shared\savedata_struct
% large 
large_dists = load('large_dists') ; large_dists = large_dists.large_dists ; 
large_inots = load('large_inots') ; large_inots = large_inots.large_inots ; 
large_latvis_inots = load('large_latvis_inots') ; large_latvis_inots = large_latvis_inots.large_latvis_inots ; 
large_mersp = load('large_mersp') ; large_mersp = large_mersp.large_mersp ; 
large_smtheeg = load('large_smtheeg') ; large_smtheeg = large_smtheeg.large_smtheeg ; 
large_sumcorrs = load('large_sumcorrs') ; large_sumcorrs = large_sumcorrs.large_sumcorrs ; 
large_peaksubs = [3,4,5,6,7,8,10,11,12,13,15,17,18,19,20,21,22,23] ; 
large_allpeaks = [64,76,74,64,66,84,72,66,74,74,70,66,64,70,64,88,70,56,84,60,66,66,54,84] ; % [21,21,20,16,16,22,16,19,16,21,17,16,16,23,15,23,21,11,19,11,18,18,19,19]*2 + 34 ;
large_postelecs = [23,56,24,57,25,58,26,59,27,60,61,62,63,64,29,30,31] ;
large_bades = [41,17,28,32,22,46] ; large_goods = zeros(1,64) ; large_goods(large_bades)=1 ; large_goods = find(large_goods==0) ; 
large_times = load('times') ; large_times = large_times.times ; large_freqs = load('freqs') ; large_freqs = large_freqs.freqs ; 
large_percs = load('large_percs'); large_percs = large_percs.large_percs; 
large_mpercs = mean(mean(large_percs(:,:,5:7),2),3); 
large_cortsz = load('large_cortsz'); large_cortsz = large_cortsz.large_cortsz; 
large_binots = load('large_binots'); large_binots = large_binots.large_binots; 

% small 
small_dists = load('small_dists') ; small_dists = small_dists.small_dists ; 
small_inots = load('small_inots') ; small_inots = small_inots.small_inots ; 
small_mersp = load('small_mersp') ; small_mersp = small_mersp.small_mersp ; 
small_smtheeg = load('small_smtheeg') ; small_smtheeg = small_smtheeg.small_smtheeg ; 
small_sumcorrs = load('small_sumcorrs') ; small_sumcorrs = small_sumcorrs.small_sumcorrs ; 
small_peaksubs = [1,2,3,4,5,6,8,9] ; 
small_allpeaks = [52,56,60,62,56,54,55,56,60] ; 
small_postelecs = [15,51,7,37,19,38,8,52,16,59,45,31,46,60,9,20,10] ; 
small_bades = [61,29,62,30,32] ; small_goods = zeros(1,64) ; small_goods(small_bades)=1 ; small_goods = find(small_goods==0) ; 
small_times = load('small_times') ; small_times = small_times.small_times ; small_freqs = load('small_freqs') ; small_freqs = small_freqs.small_freqs ; 
small_percs = load('small_percs'); small_percs = small_percs.small_percs; 
small_mpercs = mean(mean(small_percs(:,:,18:20),2),3); 
small_cortsz = load('small_cortsz'); small_cortsz = small_cortsz.small_cortsz; 

inotinds = 15:18; 
vals = load('vals') ; vals = vals.vals ; labnames = load('labnames') ; labnames = labnames.labnames ;
small_vals = load('small_vals') ; small_vals = small_vals.small_vals ; 
large_meangamma = squeeze(mean(mean(large_smtheeg(:,large_postelecs,40:90),2),3)) ;
small_meangamma = squeeze(mean(mean(small_smtheeg(:,small_postelecs,40:90),2),3)) ;
for i=1:size(vals,2)
    for j=1:size(vals,3)
        [large_valcorrs(i,j),large_valps(i,j)] = corr(vals(:,i,j),large_meangamma) ; 
        [small_valcorrs(i,j),small_valps(i,j)] = corr(small_vals(:,i,j),small_meangamma) ; 
        [large_peakcorrs(i,j),large_peakps(i,j)] = corr(vals(large_peaksubs,i,j),large_allpeaks(large_peaksubs)') ; 
        [small_peakcorrs(i,j),small_peakps(i,j)] = corr(small_vals(small_peaksubs,i,j),small_allpeaks(small_peaksubs)') ; 
        [large_peakboldcorrs(i,j),large_peakboldps(i,j)] = corr(vals(large_peaksubs,i,j),large_sumcorrs(large_peaksubs,3)) ; 
        [small_peakboldcorrs(i,j),small_peakboldps(i,j)] = corr(small_vals(small_peaksubs,i,j),small_sumcorrs(small_peaksubs,3)) ; 
        [large_peakperccorrs(i,j),large_peakpercps(i,j)] = corr(vals(large_peaksubs,i,j),large_mpercs(large_peaksubs)) ; 
        [small_peakperccorrs(i,j),small_peakpercps(i,j)] = corr(small_vals(small_peaksubs,i,j),small_mpercs(small_peaksubs)) ; 
    end
end
for i=1:16 ; disp(labnames{i,10}); end

hzs = 40:90 ; 
% distance vs power
for i=1:128 ; [large_distcorrs(i),large_distps(i)] = corr(mean(large_dists(:,large_postelecs),2),mean(mean(large_smtheeg(:,large_postelecs,i),2),3)) ; end 
for i=1:128 ; [small_distcorrs(i),small_distps(i)] = corr(mean(small_dists(:,small_postelecs),2),mean(mean(small_smtheeg(:,small_postelecs,i),2),3)) ; end 
% get all measures in all subjects
% large
large_amp = squeeze(mean(mean(large_smtheeg(:,large_postelecs,hzs),2),3)); 
large_peak = large_allpeaks'; 
large_pchange = large_mpercs;
large_nvox = large_sumcorrs(:,3);
large_dist = mean(large_dists(:,large_postelecs),2); 
large_i0 = squeeze(mean(large_inots(:,inotinds),2)); 
large_intrinsic = squeeze(mean(vals(:,[2,10],10),2)); 
large_extrinsic = squeeze(mean(vals(:,[5,13],10),2)) ;
large_calcsurf = squeeze(mean(vals(:,[1,9],20),2)); 
large_latoccsurf = squeeze(mean(vals(:,[1,9],10),2)); 
large_cortsz; 
labels = {'gamma amp.','gamma peak','BOLD %change','BOLD #voxels','distance','cancellation (i0)','intrinsic curv.','extrinsic curv.'}; 
large_vec = [large_amp,large_peak,large_pchange,large_nvox,large_dist,large_i0,large_intrinsic,large_extrinsic]; 
clear c p; 
for i=1:8
    for j=1:8
        if i==2 || j==2
            [c(i,j),p(i,j)] = corr(large_vec(large_peaksubs,i),large_vec(large_peaksubs,j)); 
        else 
            [c(i,j),p(i,j)] = corr(large_vec(:,i),large_vec(:,j)); 
        end
    end
end


for i=1:8 ; for j=1:8 ; if i<=j ; c(i,j) = 0 ; end ; end ; end
figure,subplot(2,2,1); imagesc(c(:,1:2),[-1,1]); hline(0.5:1:7.5,'k'); vline(0.5:1:7.5,'k');  colormap jet
set(gca,'XTick',1:8,'XTickLabel',labels,'YTick',1:8,'YTickLabel',labels) ; 
for i=1:8 
    for j=1:2 
        if i<=j ; text(j,i,'n/a') ; 
        else
            if p(i,j) > 0.05
                h= text(j-.25,i-.25,['rho=',num2str(round(c(i,j)*100)/100)]); 
                h = text(j-.25,i+.25,['p=',num2str(round(p(i,j)*1000)/1000)]); 
            else
                h= text(j-.25,i-.25,['* rho=',num2str(round(c(i,j)*100)/100)]); 
                h = text(j-.25,i+.25,['p=',num2str(round(p(i,j)*1000)/1000)]); 
            end
        end ; 
    end ; 
end
lcorrs = c(find(c~=0)); cl = c; 

% small 
small_amp = squeeze(mean(mean(small_smtheeg(:,small_postelecs,hzs),2),3)); 
small_peak = small_allpeaks'; 
small_pchange = small_mpercs;
small_nvox = small_sumcorrs(:,3);
small_dist = mean(small_dists(:,small_postelecs),2); 
small_i0 = squeeze(mean(small_inots(:,inotinds),2)); 
small_intrinsic = squeeze(mean(small_vals(:,[2,10],10),2)); 
small_extrinsic = squeeze(mean(small_vals(:,[5,13],10),2)) ;
small_calcsurf = squeeze(mean(small_vals(:,[1,9],20),2)); 
small_latoccsurf = squeeze(mean(small_vals(:,[1,9],10),2)); 
small_cortsz; 
labels = {'gamma amp.','gamma peak','BOLD %change','BOLD #voxels','distance','cancellation (i0)','intrinsic curv.','extrinsic curv.'}; 
small_vec = [small_amp,small_peak,small_pchange,small_nvox,small_dist,small_i0,small_intrinsic,small_extrinsic]; 
%for i=1:size(small_sumcorrs,2) ; [small_peakcorrs(i),small_peakps(i)] = corr(small_sumcorrs(small_peaksubs,i),small_allpeaks(small_peaksubs)') ; end ; 

clear c p; 
for i=1:8
    for j=1:8
        if i==2 || j==2
            [c(i,j),p(i,j)] = corr(small_vec(small_peaksubs,i),small_vec(small_peaksubs,j)); 
        else 
            [c(i,j),p(i,j)] = corr(small_vec(:,i),small_vec(:,j)); 
        end
    end
end

for i=1:8 ; for j=1:8 ; if i<=j ; c(i,j) = 0 ; end ; end ; end
subplot(2,2,2); imagesc(c(:,1:2),[-1,1]); hline(0.5:1:7.5,'k'); vline(0.5:1:7.5,'k');  colormap jet
set(gca,'XTick',1:8,'XTickLabel',labels,'YTick',1:8,'YTickLabel',labels) ; 
for i=1:8 
    for j=1:2 
        if i<=j ; text(j,i,'n/a') ; 
        else
            if p(i,j) > 0.05
                h = text(j-.25,i-.25,['rho=',num2str(round(c(i,j)*100)/100)]); 
                h = text(j-.25,i+.25,['p=',num2str(round(p(i,j)*1000)/1000)]); 
            else
                h= text(j-.25,i-.25,['* rho=',num2str(round(c(i,j)*100)/100)]); 
                h = text(j-.25,i+.25,['p=',num2str(round(p(i,j)*1000)/1000)]); 
            end
        end ; 
    end ; 
end
scorrs = c(find(c~=0)); cs = c;
[c,p] = corr(scorrs,lcorrs) ; 
subplot(2,2,3),plot(lcorrs,scorrs,'ok','LineWidth',1); lsline; 
title(['rho=',num2str(c),' p=',num2str(round(p*1000)/1000)]); xlabel('group A (rho)'); ylabel('group B (rho)'); xlim([-.7,.7]);ylim([-.85,.85]);

% for distance
figure,
subplot(2,2,1) ; 
x = large_dist ; y = large_amp ; hold on ; 
plot(x,y,'.k','MarkerSize',25) ; [rho_large,p_large] = corr(x,y) ; 
p = polyfit(x,y,1) ; r = p(1) .* x + p(2) ; plot(x,r,'-k','LineWidth',1) ; 
x = small_dist ; y = small_amp ; 
plot(x,y,'.g','MarkerSize',30) ; [rho_small,p_small] = corr(x,y) ; xlabel('distance(mm)') ; ylabel('gamma amplitude (task-rest)') ; 
p = polyfit(x,y,1) ; r = p(1) .* x + p(2) ; plot(x,r,'-g','LineWidth',1) ; 
subplot(2,2,2) ; plot(1:5,'k','LineWidth',1) ; hold on ; plot(2:6,'g','LineWidth',1) ; 
legend({['group A(n=24), ','rho=',num2str(round(rho_large*100)/100),', p=',num2str(round(p_large*1000)/1000)],...
    ['group B(n=9), rho=',num2str(round(rho_small*100)/100),', p=',num2str(round(p_small*1000)/1000)]}) ; 
subplot(2,2,3) ; plot(1:128,large_distcorrs,'k','LineWidth',1) ; hold on ; plot(1:128,small_distcorrs,'g','LineWidth',1) ; hline(0,'r') ; 
xlabel('frequency(hz)') ; ylabel('correlation (rho)') ; xlim([0,128]) ; ylim([-1,1]) ; title('*p<0.1, uncorrected') ;
for i=1:length(large_distps) ; if large_distps(i)<0.1 ; text(i,0.9,'*','Color','k') ; end ; end 
for i=1:length(small_distps) ; if small_distps(i)<0.1 ; text(i,0.85,'*','Color','g') ; end ; end 
subplot(2,2,4) ; plot(1:5,'k','LineWidth',1) ; hold on ; plot(2:6,'g','LineWidth',1) ; 
legend({['group A(n=24)'],['group B(n=9)']}) ; 

% inot vs power
figure,
for i=1:30 ; for j=1:128 ; [large_inotcorrs(i,j),large_inotps(i,j)] = corr(large_inots(:,i),mean(mean(large_smtheeg(:,large_postelecs,j),2),3)) ; end ; end
for i=1:30 ; for j=1:128 ; [small_inotcorrs(i,j),small_inotps(i,j)] = corr(small_inots(:,i),mean(mean(small_smtheeg(:,small_postelecs,j),2),3)) ; end ; end
subplot(2,2,1) ; inotinds = 15:18 ; 
x = large_i0 ; y = large_amp ; hold on ; 
plot(x,y,'.k','MarkerSize',25) ; [rho_large,p_large] = corr(x,y) ; 
p = polyfit(x,y,1) ; r = p(1) .* x + p(2) ; plot(x,r,'-k','LineWidth',1) ; 
x = small_i0 ; y = small_amp ; 
plot(x,y,'.g','MarkerSize',30) ; [rho_small,p_small] = corr(x,y) ; xlabel('I0 (A.U)') ; ylabel('gamma amplitude (task-rest)') ; 
p = polyfit(x,y,1) ; r = p(1) .* x + p(2) ; plot(x,r,'-g','LineWidth',1) ; 
subplot(2,2,2) ; plot(1:5,'k','LineWidth',1) ; hold on ; plot(2:6,'g','LineWidth',1) ; 
legend({['group A(n=24), ','rho=',num2str(round(rho_large*100)/100),', p=',num2str(round(p_large*1000)/1000)],...
    ['group B(n=9), rho=',num2str(round(rho_small*100)/100),', p=',num2str(round(p_small*1000)/1000)]}) ; 
subplot(2,2,3) ; plot(1:128,mean(large_inotcorrs(inotinds,:),1),'k','LineWidth',1) ; hold on ; 
plot(1:128,mean(small_inotcorrs(inotinds,:),1),'g','LineWidth',1) ; hline(0,'r') ; 
xlabel('frequency(hz)') ; ylabel('correlation (rho)') ; xlim([0,128]) ; ylim([-1,1]) ; title('*p<0.1, uncorrected') ;
for i=1:length(large_distps) ; if large_inotps(i)<0.1 ; text(i,0.9,'*','Color','k') ; end ; end 
for i=1:length(small_distps) ; if small_inotps(i)<0.1 ; text(i,0.85,'*','Color','g') ; end ; end 
subplot(2,2,4) ; plot(1:5,'k','LineWidth',1) ; hold on ; plot(2:6,'g','LineWidth',1) ; 
legend({['group A(n=24)'],['group B(n=9)']}) ; 

% for peak frequency
for i=1:size(large_sumcorrs,2) ; [large_peakcorrs(i),large_peakps(i)] = corr(large_sumcorrs(large_peaksubs,i),large_allpeaks(large_peaksubs)') ; end ; 
for i=1:size(small_sumcorrs,2) ; [small_peakcorrs(i),small_peakps(i)] = corr(small_sumcorrs(small_peaksubs,i),small_allpeaks(small_peaksubs)') ; end ; 
figure,subplot(2,2,1) ; 
x = large_sumcorrs(large_peaksubs,3) ; y = large_allpeaks(large_peaksubs) ; 
plot(x,y,'.k','MarkerSize',25) ; [rho_large,p_large] = corr(x,y') ; hold on ; 
p = polyfit(x,y',1) ; r = p(1) .* x + p(2) ; plot(x,r,'-k','LineWidth',1) ; 
x = small_sumcorrs(small_peaksubs,3) ; y = small_allpeaks(small_peaksubs) ; 
plot(x,y,'.g','MarkerSize',30) ; [rho_small,p_small] = corr(x,y') ; hold on ; 
p = polyfit(x,y',1) ; r = p(1) .* x + p(2) ; plot(x,r,'-g','LineWidth',1) ; ylim([45,90])
xlabel('#voxels surviving threshold') ; ylabel('peak frequency(hz)') ; 
subplot(2,2,2) ; plot(1:5,'k','LineWidth',1) ; hold on ; plot(2:6,'g','LineWidth',1) ; 
legend({['group A(n=18, removed 6 non-responders), ','rho=',num2str(round(rho_large*100)/100),', p=',num2str(round(p_large*1000)/1000)],...
    ['group B(n=8, removed 1 non-responder), rho=',num2str(round(rho_small*100)/100),', p=',num2str(round(p_small*1000)/1000)]}) ; 
% with all subjects
x = large_sumcorrs(:,3) ; y = large_allpeaks(:) ; [rho_large_all,p_large_all] = corr(x,y) ; 

% plot intrinsic curvature 
figure,
subplot(2,2,1) ; 
x = squeeze(mean(vals(:,[2,10],10),2)) ; y = squeeze(mean(mean(large_smtheeg(:,large_postelecs,hzs),2),3)) ; hold on ; 
plot(x,y,'.k','MarkerSize',25) ; [rho_large,p_large] = corr(x,y) ; 
p = polyfit(x,y,1) ; r = p(1) .* x + p(2) ; plot(x,r,'-k','LineWidth',1) ; 
x = squeeze(mean(small_vals(:,[2,10],10),2)) ; y = squeeze(mean(mean(small_smtheeg(:,small_postelecs,hzs),2),3)) ; hold on ; 
plot(x,y,'.g','MarkerSize',30) ; [rho_small,p_small] = corr(x,y) ; xlabel('intrinsic curvature (A.U)') ; ylabel('gamma amplitude (task-rest)') ; 
p = polyfit(x,y,1) ; r = p(1) .* x + p(2) ; plot(x,r,'-g','LineWidth',1) ;% xlim([5,13]) ; ylim([-.8,2.2]) ;
subplot(2,2,2) ; plot(1:5,'k','LineWidth',1) ; hold on ; plot(2:6,'g','LineWidth',1) ; 
legend({['group A(n=24), ','rho=',num2str(round(rho_large*100)/100),', p=',num2str(round(p_large*1000)/1000)],...
    ['group B(n=9), rho=',num2str(round(rho_small*100)/100),', p=',num2str(round(p_small*1000)/1000)]}) ; 

% plot mean curvature 
figure,
subplot(2,2,1) ; 
x = squeeze(mean(vals(:,[5,13],10),2)) ; y = squeeze(mean(mean(large_smtheeg(:,large_postelecs,hzs),2),3)) ; hold on ; 
plot(x,y,'.k','MarkerSize',25) ; [rho_large,p_large] = corr(x,y) ; 
p = polyfit(x,y,1) ; r = p(1) .* x + p(2) ; plot(x,r,'-k','LineWidth',1) ; 
x = squeeze(mean(small_vals(:,[5,13],10),2)) ; y = squeeze(mean(mean(small_smtheeg(:,small_postelecs,hzs),2),3)) ; hold on ; 
plot(x,y,'.g','MarkerSize',30) ; [rho_small,p_small] = corr(x,y) ; xlabel('extrinsic curvature (A.U)') ; ylabel('gamma amplitude (task-rest)') ; 
p = polyfit(x,y,1) ; r = p(1) .* x + p(2) ; plot(x,r,'-g','LineWidth',1) ;% xlim([5,13]) ; ylim([-.8,2.2]) ;
subplot(2,2,2) ; plot(1:5,'k','LineWidth',1) ; hold on ; plot(2:6,'g','LineWidth',1) ; 
legend({['group A(n=24), ','rho=',num2str(round(rho_large*100)/100),', p=',num2str(round(p_large*1000)/1000)],...
    ['group B(n=9), rho=',num2str(round(rho_small*100)/100),', p=',num2str(round(p_small*1000)/1000)]}) ; 

figure,
m_largemersp = squeeze(mean(mean(large_mersp(:,:,20:45,50:170),3),4)) ; 
m_smallmersp = squeeze(mean(mean(small_mersp(:,:,20:45,50:170),3),4)) ; 
[c,p] = corr(m_largemersp) ; 
for i=1:size(c,1);for j=1:size(c,2); if i<=j ; c(i,j) = 0; end; end ;end
for i=1:size(p,1);for j=1:size(c,2); if i<=j ; p(i,j) = 0.05; end; end ;end
subplot(2,2,3); imagesc(c,[-1,1]);
subplot(2,2,4); imagesc(p,[0,0.1]);
subplot(2,2,1) ; 
x = m_largemersp(:,1) ; y = m_largemersp(:,5) ; 
plot(x,y,'.k','MarkerSize',25) ; [rho_large,p_large] = corr(x,y) ; hold on ; 
p = polyfit(x,y,1) ; r = p(1) .* x + p(2) ; plot(x,r,'-k','LineWidth',1) ; 
x = m_smallmersp(:,1) ; y = m_smallmersp(:,2) ; 
plot(x,y,'.g','MarkerSize',30) ; [rho_small,p_small] = corr(x,y) ; hold on ; xlabel('stimulus 1 (unperturbed)') ; ylabel('stimulus 5 (10% randomized)') ;
p = polyfit(x,y,1) ; r = p(1) .* x + p(2) ; plot(x,r,'-g','LineWidth',1) ; 
subplot(2,2,2) ; plot(1:5,'k','LineWidth',1) ; hold on ; plot(2:6,'g','LineWidth',1) ; 
legend({['group A(n=24), ','rho=',num2str(round(rho_large*100)/100),', p=',num2str(round(p_large*1000)/1000)],...
    ['group B(n=9), rho=',num2str(round(rho_small*100)/100),', p=',num2str(round(p_small*1000)/1000)]}) ; 

figure,
for i=1:6 ; 
    subplot(2,6,i) ; imagesc(large_times,large_freqs,squeeze(mean(large_mersp(:,i,:,:),1)),[-2,2]) ; axis xy ; 
    if i==1 ; xlabel('time(sec)') ; ylabel('frequency(hz)') ; end
end
for i=1:3 ; 
    subplot(2,6,i+6) ; imagesc(small_times,small_freqs,squeeze(mean(small_mersp(:,i,:,:),1)),[-3,3]) ; axis xy ; 
    if i==1 ; xlabel('time(sec)') ; ylabel('frequency(hz)') ; end
end
subplot(2,6,10) ; imagesc(large_times,large_freqs,squeeze(mean(mean(large_mersp,1),2)),[-1,1]); axis xy;  hline([40,90],'k'); vline([0,2],'k'); 



figure,
subgamma = (squeeze(mean(mean(mean(large_mersp(:,[1,5],25:40,80:150),2),3),4))) ;  
[~,si] = sort(subgamma,'descend') ; 
for i=1:length(si) ; 
    subplot(2,12,i) ; imagesc(large_times,large_freqs,squeeze(mean(large_mersp(si(i),[1],:,1:end-10),2)),[-3,3]) ; axis xy ; 
    title(['subject ',num2str(i)]); if i==1 ; xlabel('time(s)') ;ylabel('frequency(hz)'); end
    %text(45,55,['subject ',num2str(i)]) ;
    %set(gca,'XTick',[],'YTick',[]) ; 
end
figure,
for i=1:6 ; subplot(2,6,i) ; imagesc(large_times,large_freqs,squeeze(large_mersp(si(2),i,:,:)),[-3,3]) ; axis xy; end
for i=1:6 ; subplot(2,6,i+6) ; imagesc(large_times,large_freqs,squeeze(large_mersp(si(22),i,:,:)),[-3,3]) ; axis xy; end
figure,
subgamma = (squeeze(mean(mean(mean(small_mersp(:,[1,2],20:40,80:150),2),3),4))) ;  
[~,si] = sort(subgamma,'descend') ; 
for i=1:length(si)
    subplot(2,12,i) ; imagesc(small_times,small_freqs,squeeze(mean(small_mersp(si(i),[1],:,1:end-10),2)),[-4,4]) ; axis xy ; 
    title(['subject ',num2str(i)]) ; if i==1; xlabel('time(s)'); ylabel('frequency(hz)'); end
end
figure,
subplot(2,3,1) ; imagesc([-1,1]); colorbar; subplot(2,3,2); imagesc([-2,2]) ; colorbar ; subplot(2,3,3); imagesc([-3,3]) ; colorbar; 
subplot(2,3,4) ; imagesc([-4,4]); colorbar; 

figure,
plot(squeeze(mean(mean(large_mersp(6,[1,5],15:end,80:150),2),4)),'b','LineWidth',1) ; hold on ;
plot(squeeze(mean(mean(large_mersp(18,[1,5],15:end,80:150),2),4)),'r','LineWidth',1) ;
set(gca,'XTick',1:5:45,'XTickLabel',30:10:120) ; ylabel('amplitude (db)'); xlabel('frequency(hz)'); 
legend({'subject 9','subject 6'}); 



cd C:\Users\butr2901\Desktop\structpap
