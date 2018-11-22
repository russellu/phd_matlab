cd('C:\shared\badger\EEG-fMRI 2015-10-01\genevieve') ; 
clear all ; close all ; 
mongs = dir('*vhdr') ; 
for i=1:length(mongs) ; 
   EEG = pop_loadbv('.',mongs(i).name) ;  
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
   EEG = denoise_grad(EEG) ;  EEG = pop_resample(EEG,512) ; EEG = denoise_bcg(EEG) ;    
   eegs{i} = EEG ; 
   if i==1 ; merged = EEG  ;else merged = pop_mergeset(EEG,merged) ; end     
end
cd('C:\shared\badger\EEG-fMRI 2015-10-01\') ; regs = dir('reg_*') ; 
for r=1:length(regs) ; fmris(r) = load_untouch_nii(regs(r).name) ; end

e1 = eegs{1} ; 
volOnsets = find(strcmp('R128',{e1.urevent.type})) ; 
trignums = 1:16 ; trigs = {} ; 
for i=1:length(trignums) ; if trignums(i) < 10 ; trigs{i} = ['S  ',num2str(i)] ; else trigs{i} = ['S ',num2str(i)] ; end ; end










