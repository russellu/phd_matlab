perc = load('c:/shared/all_white_normals/fmris/sub_vincent/subcorrs') ; perc = perc.subcorrs ; 
vox = load('c:/shared/all_white_normals/fmris/sub_vincent/meanvoxcorrs') ; vox = vox.meanvoxcorrs ; 
dist = load('c:/shared/all_white_normals/a1_good/sub_vincent/corrmat.mat') ; dist = dist.corrmat ; 
curv = load('c:/shared/all_white_normals/a1_good/sub_vincent/r.mat') ; curv = curv.r ; 

plot(curv(20,:),'r') ; hold on ; 
plot(dist(:,55),'k') ; hline(0,'k') ; 
plot(vox,'m') ; 
plot(perc,'b') ; ylabel('correlation(rho)') ; 
xlabel('frequency(hz)') ; set(gca,'XTick',1:5:60,'XTickLabel',1:10:120) ; 
legend({'curvature','distance','#voxels','%change'}) ; 