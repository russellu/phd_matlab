% estimate hrf from event related designclear all ; close all ; 
% stimulus starts in the bottom left quadrant (225deg=S1). 
trigs = {'S  1','S  2','S  3'} ; 
subs= {'alex','dina','genevieve','jeremie','karl','russell','valerie'} ;

for sub=1:length(subs)
sub = subs{sub} ; 
cd(['c:/shared/badger_mri/',sub,'/nii']) ; ls 
fmris=dir('reg_mc*gamma*nii*') ; 
cd(['c:/shared/badger_eeg/',sub]) ; 
eegs = dir('1hz*gamma*set') ; 

stimcounts = [1,1,1] ;  
allstims = zeros(size(fimg,1),size(fimg,2),size(fimg,3),round(15./0.693)) ; 
for f=1:length(fmris) ; 
   cd(['c:/shared/badger_mri/',sub,'/nii']) ; ls 
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
   trtopeak = find(spm_hrf(TR)==max(spm_hrf(TR))) ; 
   for i=1:length(trigs) ; % for all orientations
          trigi = find(strcmp(trigs{i},events)) ; 
          for j=1:length(trigi)
              latij = cell2mat(lats(trigi(j))) ; 
              diffs = onsetlats-latij ; 
              indij = find(abs(diffs) == min(abs(diffs))) ; 
              allstims = allstims + fimg(:,:,:,indij-round(5/TR):indij+round(10/TR)) ; 
          end
   end
end
cd(['c:/shared/badger_mri/',sub,'/nii']) ; ls 
ref = load_untouch_nii('hrf_ref.nii.gz') ; 
ref.img = allstims ; save_untouch_nii(ref,'allstims.nii.gz') ; 

end
