clear all ; close all ; 
cd c:/shared/lead ; 
lead = load('lead') ; lead = lead.lead ; 
locs = lead.GridLoc ; 
%Vq = griddata(X,Y,Z,V,Xq,Yq,Zq)

[svx,six] = sort(locs(:,1),'descend') ; ux = unique(svx) ;
[svy,siy] = sort(locs(:,2),'descend') ; uy = unique(svy) ; 
[svz,siz] = sort(locs(:,3),'descend') ; uz = unique(svz) ; 

for i=1:size(locs,1)
   ix = find(locs(i,1)==ux) ; 
   iy = find(locs(i,2)==uy) ; 
   iz = find(locs(i,3)==uz) ; 
   grid(ix,iy,iz) = 1 ; 
end

wcell = load('wcell')  ; wcell = wcell.wcell ; 
winv = pinv(wcell{1}*wcell{2}) ; 
cd(['C:\shared\simdenoise\alex']) ;
EEG = pop_loadbv('.','retino_gamma_01.vhdr') ; 
EEG = pop_chanedit(EEG,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ; 

leadfield = lead.Gain ; 
leadfield(isnan(leadfield)) = 0 ; 
[c,p] = corr(winv,leadfield) ; 
cx = c(:,1:size(c,2)/3) ; 
cy = c(:,size(c,2)/3+1:size(c,2)/3+size(c,2)/3) ; 
cz = c(:,end-size(c,2)/3+1:end) ; 
cx = c(:,1:3:end) ; cy = c(:,2:3:end) ; cz = c(:,3:3:end) ; 

cm = (cx.^2)+(cy.^2)+(cz.^2) ; 
clear compcorr
for c=1:64 ; disp(c) ; 
for i=1:size(locs,1)
   ix = find(locs(i,1)==ux) ; 
   iy = find(locs(i,2)==uy) ; 
   iz = find(locs(i,3)==uz) ; 
   %grid(ix,iy,iz) = 1 ; 
   compcorr(ix,iy,iz,c) = cm(c,i) ; 
end
end


cd c:/shared/lead ; 
save_nii(make_nii(compcorr),'compcorr.nii.gz') ; 
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(fliplr(squeeze(max(compcorr(:,:,:,i),[],3)))) ; axis xy ; title(i) ; end
figure,for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),EEG.chanlocs) ; title(i) ; end


