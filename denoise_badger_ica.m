clear all ; close all; 

subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'};
bades = {[1,4,26],[1,2,4,6,16,17,19,23],[1,2,4,6,19],[3,14],[1,4,6,],[2,4,7,19],[1,5,8,10,23,29,28],[1,4,8,10,11,13]};

for sb=1:length(subs)
cd(['E:\badger_eeg\',subs{sb}]);

merged = pop_loadset('merged_denbcg.set'); 

fullcomps = load('fullcomps'); fullcomps = fullcomps.fullcomps; weights = fullcomps{1}; sphere = fullcomps{2}; 
acts = weights*sphere*merged.data; 
winv = pinv(weights*sphere); 


end




