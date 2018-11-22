clear all  ; close all
subs = {'charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','lisa','marc','marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; ; 
for s=1:length(subs)
cd(['C:\shared\resmerged\',subs{s}]) ; 
postersp = load('postersp') ; postersp = postersp.postersp(:,1:128) ; 
bdersp = load('dersp') ; bdersp = bdersp.dersp ; 
melecs = squeeze(mean(bdersp(:,:,40:150),3)) ; allmelecs(s,:,:) = melecs ; 
%melecs = postersp ; 

if s==1 ;
EEG = pop_loadset('merged.set') ; 
EEG = pop_chanedit(EEG,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ; 
end
cd(['C:\shared\all_white_normals\a1_good\sub_',subs{s}]) ;
corrs = load_untouch_nii('meancorrs_t1.nii.gz') ;
corrs.img(:,1:180,:) = 0 ; 
%allmasks(:,:,:,s) = corrs.img;  
locs = load_untouch_nii('elecbrain.nii.gz') ; 
length(unique(locs.img(:)))
clear cx cy cz
for i=1:65
   [cx(i),cy(i),cz(i)] = centmass3(locs.img==i) ;  
    
end
[rx,ry,rz] = centmass3(corrs.img>.35) ; 
figure,imagesc(sum(corrs.img>.35,3)) ; 
order = load('c:/shared/all_white_normals/a2_good/elecorder.mat') ; order = order.elecorder ; 
labs = {EEG.chanlocs.labels} ; clear allinds 
bads = zeros(1,length(labs)) ; 
for i=1:length(labs)
    ind = find(strcmpi(labs{i},order)) ; 
    if ~isempty(ind)
    allinds(i) = ind ;  
    else 
        allinds(i) = 1  ; bads(i) = 1 ; 
    end
end ; goods = find(bads==0) ; bads = find(bads==1) ; 
allelecs(s,:,:) = melecs ; 
mdx = cx(allinds) ; mdy = cy(allinds) ; mdz = cz(allinds) ; 

allcoords(s,1,:) = mdx ; allcoords(s,2,:) = mdy ; allcoords(s,3,:) = mdz ; 

diffs = [mdx-rx;mdy-ry;mdz-rz] ;
sqrdiff = sqrt(sum(diffs.^2,1)) ; 

sqrdiff(bads) = [] ; 
melecs(bads,:) = [] ; 
alldiffs(s,:) = sqrdiff ; %allelecs(s,:,:) = melecs ; 
c = corr(sqrdiff',melecs) ;
cs(s,:) = c ; 
end

newelecs = allelecs ; newelecs(:,bads,:) = [] ; clear corrs ;

clear corrs ; 
for i=1:60
    for j=1:62
         corrs(i,j) = corr2(squeeze(alldiffs(:,j)),squeeze(newelecs(:,j,i))) ; 
    end
end
topocorrs = zeros(size(corrs,1),64) ; 
topocorrs(:,goods) = corrs(:,:) ; 


shadedErrorBar([],mean(cs,1),std(cs,0,1)/sqrt(22)) ; hline(0,'k') ; set(gca,'XTick',0:5:60,'XTickLabel',0:10:120) ; xlabel('frequency(hz)') ; ylabel('corr(rho)') ; 
for i=1:60 ; [h,p,ci,stats] = ttest(cs(:,i)) ; allts(i) = stats.tstat ; allps(i) = p ; end ; 
for i=1:60 ; if allps(i) < 0.05/60 ; text(i,.55,'*') ; end ; end ; set(gca,'XTick',0:5:60,'XTickLabel',0:10:120) ; xlabel('frequency(hz)') ; ylabel('corr(rho)') ; 


subplot(1,2,1) ; topoplot(squeeze(mean(mean(allmelecs(:,:,5:12),1),3)),EEG.chanlocs,'maplimits',[-1.5,1.5],'electrodes','off') ;
subplot(1,2,2) ; topoplot(squeeze(mean(mean(allmelecs(:,:,30:40),1),3)),EEG.chanlocs,'maplimits',[-.5,.5],'electrodes','off') ;

subplot(1,2,1) ; imagesc(g2,[-.5,.5]) ; axis xy ; subplot(1,2,2) ; imagesc(g1,[-1.5,1.5]) ; axis xy ;  






