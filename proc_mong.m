clear all; close all; 
cd c:/shared/mong_eeg/
cd MONG_01_RB; 

bvs = dir('bv*vhdr'); 
for bv=1:length(bvs)
   eeg = pop_loadbv('.',bvs(bv).name);  
   eeg = pop_chanedit(eeg,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
   if bv==1; merged = eeg; else merged = pop_mergeset(eeg,merged); end

end
filtmerged = eegfiltfft(merged.data,merged.srate,1,90); 
[weights,sphere] = runica(filtmerged,'maxsteps',128);
winv = pinv(weights*sphere); 
for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),merged.chanlocs) ; end




