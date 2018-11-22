clear all ; close all ; 
% stimulus starts in the bottom left quadrant (225deg=S1). 
trigs{1} = {'S 11','S 12','S 13','S 14'} ; 
trigs{2} = {'S 21','S 22','S 23','S 24'} ; 
trigs{3} = {'S 31','S 32','S 33','S 34'} ; 
trigs{4} = {'S 41','S 42','S 43','S 44'} ; 
trigs{5} = {'S 51','S 52','S 53','S 54'} ; 
trigs{6} = {'S 61','S 62','S 63','S 64'} ; 
subs= {'alex','dina','genevieve','jeremie','karl','russell','tegan','valerie'} ;

for ss=6:length(subs) 
sub = subs{ss} ; 
cd(['c:/shared/badger_mri/',sub,'/nii']) ; ls 
fmris=dir('bp_*allstims*nii*') ; 
cd(['c:/shared/badger_eeg/',sub]) ; 
eegs = dir('allfreq*allstim*set') ; 
% just flip the x component of the vector

startangles = [0,60,120,180,240,300] ;

if ss==1 || ss==4
    realangles = mod(225-startangles,360) ; 
    r = ones(1,6) ; 
    radangles = (realangles*pi)/180 ; 
    [x,y] = pol2cart(radangles,r) ; 
    [theta,rho] = cart2pol(-x,y) ; 
    deg = theta*180/pi ; 
    realangles = mod(deg,360) ;    
else
    realangles = mod(225-startangles,360) ; 
end

clear allbrains
for f=1:length(fmris) ; 
   cd(['c:/shared/badger_mri/',sub,'/nii']) ; ls 
   delay = load_untouch_nii('delay_mean.nii.gz') ; 
   delayimg = delay.img(:,:,:,1,1) ; 
   delayimg(delayimg>7) = 0 ; % set unrealistic delays to 0
   fmri = load_untouch_nii(fmris(f).name) ; 
   fimg = fmri.img ; 
   mask = load_untouch_nii('epi_v1_cluster.nii.gz') ; 
   maskim = imdilate(mask.img,strel(ones(13,13,13))) ; 
   maskinds = find(maskim==1) ;
   for maskind=1:length(maskinds)
       [ix,iy,iz] = ind2sub(size(maskim),maskinds(maskind)) ;  
       postvox(maskind,:) = squeeze(fimg(ix,iy,iz,:)) ;
   end
    
 

   TR = fmri.hdr.dime.pixdim(5) ; 
   cd(['c:/shared/badger_eeg/',sub]) ; ls 
   eeg = pop_loadset(eegs(f).name) ; 
   
    trigsr = {eeg.urevent.type} ; latsr = cell2mat({eeg.urevent.latency}) ; 
    r128s = find(strcmp(trigsr,'R128')) ;
    rlats = latsr(r128s) ; 
    eeg = pop_select(eeg,'nopoint',[1,rlats(1) ; rlats(length(rlats)),size(eeg.data,2)]) ; 

   
   events = {eeg.urevent.type} ; 
   onsetinds = find(strcmp('R128',events)) ; 
   lats = {eeg.urevent.latency} ; 
   onsetlats = cell2mat(lats(onsetinds)) ; 
   vlats = cell2mat(lats(onsetinds)) ; 

   angleincr = 360/round(10/TR) ; 
   clear fullangles
   for i=1:length(trigs) ; 
       if ss==1 || ss==4
           fullangles(i,:) = mod(realangles(i):angleincr:realangles(i)+390,360) + rand(1,16) ;
       else
           fullangles(i,:) = mod(realangles(i)+360:-angleincr:realangles(i),360) + rand(1,15) ;
       end
   end
   uniques = unique(fullangles) ; 
   allangles = zeros(size(fimg,1),size(fimg,2),size(fimg,3),length(uniques)) ; 
   for i=1:length(trigs) % for all starting locations
       for j=1:length(trigs{i}) ; % for all orientations
            trigij = find(strcmp(trigs{i}(j),events)) ; 
            latij = cell2mat(lats(trigij)) ; 
            diffs = onsetlats-latij ; 
            indij = find(abs(diffs) == min(abs(diffs)),1) ; 
            anglesi = fullangles(i,:) ; 
            for a=1:length(anglesi)
                uniqueind = find(uniques==anglesi(a)) ; 
               for maskind=1:length(maskinds)
                   [ix,iy,iz] = ind2sub(size(maskim),maskinds(maskind)) ;  
                   ts = squeeze(fimg(ix,iy,iz,:)) ; 
                   allangles(ix,iy,iz,uniqueind) = allangles(ix,iy,iz,uniqueind) + ...
                       mean(fimg(ix,iy,iz,indij+a+round(delayimg(ix,iy,iz))-1:indij+a+round(delayimg(ix,iy,iz))+1),4) ;% - mean(fimg(:,:,:,indij+a-5:indij+a-6),4) ;
               end            
            end
       end
   end
 
   % get the stimulus onset latency and compare it to the TR latencies.
   % convert the stimulus time into TRs. 
   % how to find things oyu want: convolve the design with the HRF
   allbrains(f,:,:,:,:) = allangles ; 
end

cd(['c:/shared/badger_mri/',sub,'/nii']) ; ls  
angleref = load_untouch_nii('angleref_90.nii.gz') ; 

mbrains = squeeze(mean(allbrains,1)) ; 
smoothbrains = zeros(size(mbrains)) ; 

for i=1:size(mbrains,1)
    for j=1:size(mbrains,2)
        smoothbrains(i,j,:,:) = imfilter(squeeze(mbrains(i,j,:,:)),fspecial('gaussian',[1,5],3)) ; 
    end
end

angleref.img = squeeze(smoothbrains) ; 
save_untouch_nii(angleref,'allangles.nii.gz') ;

anat = load_untouch_nii('fref.nii.gz') ; 
zinds = zeros(size(smoothbrains,1),size(smoothbrains,2),size(smoothbrains,3)) ; 
for i=1:size(smoothbrains,1)
    for j=1:size(smoothbrains,2)
        for k=1:size(smoothbrains,3)
            zinds(i,j,k) = find(squeeze(smoothbrains(i,j,k,:))==max(squeeze(smoothbrains(i,j,k,:))),1) ;             
        end
    end
end
anat.img = zinds ; 
save_untouch_nii(anat,'angleinds.nii.gz') ; 

end














