clear all ; close all ; 
% stimulus starts in the bottom left quadrant (225deg=S1). 
trigs{1} = {'S 11','S 12','S 13','S 14'} ; 
trigs{2} = {'S 21','S 22','S 23','S 24'} ; 
trigs{3} = {'S 31','S 32','S 33','S 34'} ; 
trigs{4} = {'S 41','S 42','S 43','S 44'} ; 
trigs{5} = {'S 51','S 52','S 53','S 54'} ; 
trigs{6} = {'S 61','S 62','S 63','S 64'} ; 

subs= {'alex','dina','genevieve','jeremie','karl','russell','valerie'} ;
for ss=4%:length(subs) 
sub = subs{ss} ; 
peaks_TR = [7,7,7,7,7,7,7] ;
cd(['c:/shared/badger_mri/',sub,'/nii']) ; ls 
fmris=dir('blur_bp_reg_mc*allstims*nii*') ; 
cd(['c:/shared/badger_eeg/',sub]) ; 
eegs = dir('1hz*allstim*set') ; 

startangles = 180-[0,60,120,180,240,300] ;
realangles = mod(225-startangles,360) ; 


for f=1:length(fmris) ; 
   cd(['c:/shared/badger_mri/',sub,'/nii']) ; ls 
   delay = load_untouch_nii('delay_avg.nii.gz') ; 
   delayimg = delay.img(:,:,:,1,1) ; 
   
   
   fmri = load_untouch_nii(fmris(f).name) ; 
   fimg = fmri.img ; 

   TR = fmri.hdr.dime.pixdim(5) ; 
   cd(['c:/shared/badger_eeg/',sub]) ; ls 
   eeg = pop_loadset(eegs(f).name) ; 
   events = {eeg.urevent.type} ; 
   onsetinds = find(strcmp('R128',events)) ; 
   lats = {eeg.urevent.latency} ; 
   onsetlats = cell2mat(lats(onsetinds)) ; 
   vlats = cell2mat(lats(onsetinds)) ; 
   trtopeak = peaks_TR(f) ; 
   %boldtrs = round(1:18/TR) ; 
   angleincr = 360/round(10/TR) ; 
   for i=1:length(trigs) ; 
       %fullangles(i,:) = mod(realangles(i)+360:-angleincr:realangles(i),360) + rand(1,15) ;
       fullangles(i,:) = mod(realangles(i):angleincr:realangles(i)+390,360) + rand(1,16) ;
   end
   uniques = unique(fullangles) ; 
   allangles = zeros(size(fimg,1),size(fimg,2),size(fimg,3),length(uniques)) ; 
   for i=1:length(trigs) % for all starting locations
       for j=1:length(trigs{i}) ; % for all orientations
            trigij = find(strcmp(trigs{i}(j),events)) ; 
            latij = cell2mat(lats(trigij)) ; 
            diffs = onsetlats-latij ; 
            indij = find(abs(diffs) == min(abs(diffs))) ; 
            anglesi = fullangles(i,:) ; 
            for a=1:length(anglesi)
                uniqueind = find(uniques==anglesi(a)) ; 
                allangles(:,:,:,uniqueind) = allangles(:,:,:,uniqueind) + mean(fimg(:,:,:,indij+a+trtopeak-1:indij+a+trtopeak+1),4) ;% - mean(fimg(:,:,:,indij+a-5:indij+a-6),4) ;
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

anat = load_untouch_nii('anat.nii.gz') ; 
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














