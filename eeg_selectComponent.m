clear all ; close all
cd c:/shared/papsaves ; ls 
subs = {'sub_alex','sub_charest','sub_esteban','sub_fabio','sub_gab','sub_gabriella','sub_genevieve','sub_gina','sub_jeremie','sub_julie','sub_katrine','sub_lisa'...
        ,'sub_marc','sub_marie','sub_mathieu','sub_maxime','sub_mingham','sub_patricia','sub_po','sub_russell','sub_sunachakan','sub_vincent'} ;
stimnames = {'5%cont','33%cont','unpt','10%rnd','60%rnd'} ; stims = [2,3,1,5,6] ; 
EEG = pop_loadset('ica_notch85.set') ; 
elecs = load('elecact.mat') ; elecs = elecs.elecact ; 
comps = load('icaact.mat') ; comps = comps.icaact ; 
winvs = load('allwinvs.mat') ; winvs = winvs.allwinvs ; 
freqs = load('freqs.mat') ; freqs = freqs.freqs ; 
times = load('times.mat') ; times = times.times ;
wmask = double(squeeze(mean(mean(mean(mean(elecs(:,:,:,(freqs<27 & freqs>15),(times<2000&times>0)),1),2),4),5))).^2 ; 
t_roi = find(times>500 & times<1500) ; 
f_roi = find(freqs>40 & freqs<100) ; 
post = [23,56,24,57,25,58,26,59,27,60,61,62,63,64,29,30,31] ; 
elecroi = squeeze(mean(mean(elecs(:,:,post,f_roi,t_roi),1),3)) ; 
comprois = squeeze(comps(:,:,:,f_roi,t_roi)) ; 
for i=1:22 ;
    for j=1:6
        for k=1:64
            compcorrs(i,j,k) = corr2(squeeze(elecroi(j,:,:)),squeeze(comprois(i,j,k,:,:))) ;         
        end
    end
end
meancompcorrs = squeeze(mean(compcorrs,2)) ; 

f_roi = find(freqs>5 & freqs<40) ; 
post = [23,56,24,57,25,58,26,59,27,60,61,62,63,64,29,30,31] ; 
elecroi = squeeze(mean(mean(elecs(:,:,post,f_roi,t_roi),1),3)) ; 
comprois = squeeze(comps(:,:,:,f_roi,t_roi)) ; 
for i=1:22 ;
    for j=1:6
        for k=1:64
            compcorrs2(i,j,k) = corr2(squeeze(elecroi(j,:,:)),squeeze(comprois(i,j,k,:,:))) ;         
        end
    end
end
meancompcorrs2 = squeeze(mean(compcorrs2,2)) ; 

for i=1:size(winvs,1)
    for j=1:size(winvs,3)
        winvcorrs(i,j) = corr2(wmask,squeeze(abs(winvs(i,:,j)))') ; 
    end
end
finalcorrs = (meancompcorrs*2.5+meancompcorrs2/2)+((winvcorrs>0).*winvcorrs.^2) ; 
[sv,si] = sort(finalcorrs,2,'descend') ; 
for i=1:22 ; figure ; for j=1:30 ; subplot(5,6,j) ; imagesc(squeeze(comps(i,1,si(i,j),:,:)),[-3,3]) ; title(winvcorrs(i,si(i,j))); end ; suptitle(subs{i}) ;end 
clear sorties
for i=1:22 ; sorties(i,:,:,:) = squeeze(mean(comps(i,:,si(i,1:2),:,:),3)) ; end
msorts = squeeze(mean(mean(mean(sorties(:,[1,5],freqs<80 & freqs>50,times>0&times<2000),2),3),4)) ; 
[sv,sig] = sort(msorts,'descend') ; 
stims = [1,3,2,5,6] ; 

for i=1:22 ; 
    h = subplot(5,22,i) ; 
    imagesc(flipud(squeeze((sorties(sig(i),1),:,:)))),[-3,3]) ; set(gca,'XTick',[],'YTick',[]) ; 
     p = get(h, 'pos') ; 
    if i>1 
        p(1) = prevleft + prevwidth ; 
        set(h, 'pos', p);
    end
    prevleft = p(1) ; 
    prevwidth = p(3) ; 
end


for i=1:22
    for j=1:5 
        ssorts(i,j,:,:) = flipud(squeeze(ssorts(i,j,:,:))) ; 
        
    end
end

subs = {'sub_alex','sub_charest','sub_esteban','sub_fabio','sub_gab','sub_gabriella','sub_genevieve','sub_gina','sub_jeremie','sub_julie','sub_katrine','sub_lisa'...
        ,'sub_marc','sub_marie','sub_mathieu','sub_maxime','sub_mingham','sub_patricia','sub_po','sub_russell','sub_sunachakan','sub_vincent'} ;

cd c:/shared/regf1; ls 
for s=1:length(subs)
   a = load_untouch_nii(['reg_',subs{s},'.nii.gz']) ; 
   anats(s,:,:,:) =  a.img ;
   c = load_untouch_nii(['comps_',subs{s},'.nii.gz']) ; 
   comps(s,:,:,:) = c.img ; 
end
for i=1:22 ; 
    h = subplot(5,22,i) ; 
    imagesc(flipud(squeeze((sorties(sig(i),1),:,:)))),[-3,3]) ; set(gca,'XTick',[],'YTick',[]) ; 
     p = get(h, 'pos') ; 
    if i>1 
        p(1) = prevleft + prevwidth ; 
        set(h, 'pos', p);
    end
    prevleft = p(1) ; 
    prevwidth = p(3) ; 
end










