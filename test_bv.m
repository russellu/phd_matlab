%{
cd c:/shared/eeg_tests 
ls
EEG = pop_loadset('retino_gamma_01_Pulse Artifact Correction.set') ; 
filt = eegfiltfft(EEG.data,EEG.srate,2,250) ; EEG2 = EEG ; EEG2.data = filt ; 
ica = pop_runica(EEG2,'runica') ; 
%}
clear all ; close all
cd c:/shared/badger_eeg/tegan ; 
grads = dir('subgrad*set') ; 
for i=1:length(grads)
    EEGi = pop_loadset(grads(i).name) ;  
    [datazs{i},eeg4s{i}] = fbcg(EEGi) ; 
end
l = 0 ; 
for i=1:length(datazs)
    l = l + size(datazs{i},2) ; 
end
alldatazs = zeros(64,l) ; 
l=1 ; 
for i=1:length(datazs)
    alldatazs(:,l:l+size(datazs{i},2)-1) = datazs{i} ; 
    l = l+size(datazs{i},2) ; 
end

EEGi = pop_chanedit(EEGi,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
%topoplot(alldatazs(:,ind),EEGi.chanlocs,'electrodes','numbers') ; 

motions = mean(alldatazs([62,12,54,56],:)) - mean(alldatazs([61,29,55,13],:)) ; 
smoothmotions = smooth(motions.^2,100) ; 
[sv,si] = sort(smoothmotions,'descend') ; 
bads = si(1:round(length(si)/8)) ; 

newdatazs = alldatazs ; newdatazs(:,bads) = [] ; 
[weights,sphere] = runica(newdatazs,'maxsteps',128) ; 
winv = pinv(weights*sphere) ; acts = (weights*sphere)*newdatazs ; 
for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),EEGi.chanlocs) ; end 
[s,f] = spectopo(acts,0,250,'plot','off') ; 

for i=1:length(eeg4s)
    pop_saveset(eeg4s{i},['bcgmot_',grads(i).name]) ;
end
icadat = weights*sphere ; 
save('icadat','icadat')  ;


%{
cd c:/shared/badger_eeg/tegan/outside ; ls
out = dir('*outside*vhdr') ; out = pop_loadbv('.',out.name) ; 
out = pop_resample(out,250) ;
out.data = eegfiltfft(out.data,250,1,250) ; 
outica = pop_runica(out,'runica') ; 
figure,for i=1:64 ; subplot(5,13,i) ; topoplot(outica.icawinv(:,i),EEGi.chanlocs) ; end
%}

for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),EEGi.chanlocs) ; title(i) ; end 

eeg = eeg4s{6} ; 

eact = icaact(eeg.data,weights*sphere,0) ; 




