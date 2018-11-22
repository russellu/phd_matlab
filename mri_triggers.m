cd('c:/Vision/Raw Files') ; ls 
clear all ; close all ; 
ssveps = dir('Test_Russell_2015-08-13*vhdr') ; 

for f=1:length(ssveps) ; 
    EEG = pop_loadbv('.',ssveps(f).name) ; 
    EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
    EEG = pop_resample(EEG,500) ; 
    eegs{f} = EEG ; 
end

%%%% getting the MR triggers: assuming the FMRI time points are a subset of
%%%% the EEG time-points, you can just extract the FMRI epochs at which the 
fvols = find(strcmp('R128',{eegs{1}.urevent.type})) ; 
lats = cell2mat({eegs{1}.urevent.latency}) ; 
fvol_lats = lats(fvols) ; 

trigs = {'S  1','S  2','S  3','S  4'} ; 
for t=1:length(trigs)
   triglats_t = find(strcmp(trigs{t},{eegs{1}.urevent.type})) ; 
   lats_t = lats(triglats_t) ; 
   for i=1:length(lats_t)
       latdiffs = abs(fvol_lats - lats_t(i)) ; 
       trigvols(t,i) = find(latdiffs == min(latdiffs)) ; 
   end
end
cd c:/shared/UTE ; ls ; fcount = 1 ; clear allepochs; 
for f=4:6
fmri = load_untouch_nii(['MONG_01_RB_WIP_fMRI_S2.1_TR115_SENSE_',num2str(f),'_1.nii']) ; 
fim = fmri.img ; clear fepochs  
for i=1:size(trigvols,2)
    fepochs(:,:,:,i,1:51) = fim(:,:,:,trigvols(2,i)-10:trigvols(2,i)+40) ; 
end
allepochs(fcount,:,:,:,:,:) = fepochs ; fcount = fcount + 1 ; 
end
mepochs = squeeze(mean(mean(allepochs,1),5)) ; 
save_nii(make_nii(mepochs),'mepochs.nii.gz') ; 




