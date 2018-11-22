clear all ; close all; 

subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'};
for sb=1:length(subs)
cd(['E:\badger_eeg\',subs{sb}]);

%{
denbcgs = dir('denbcg*set'); 
for st=1:length(denbcgs)
   disp(denbcgs(st).name);  
   if st==1; merged = pop_loadset(denbcgs(st).name); else merged = pop_mergeset(merged,pop_loadset(denbcgs(st).name)); end
end

pop_saveset(merged,'merged_denbcg.set'); 
%}
merged = pop_loadset('merged_denbcg.set'); 

%{
filtmerged = eegfiltfft(merged.data,merged.srate,1,128);
[weights,sphere] = runica(filtmerged,'maxsteps',128); 
winv = pinv(weights*sphere) ; 

fullcomps{1} = weights; fullcomps{2} = sphere; save('fullcomps','fullcomps'); 

figure,for i=1:64  ; subplot(5,13,i) ; topoplot(mat2gray(winv(:,i)),merged.chanlocs) ; end ;suptitle(subs{sb}); 

acts = weights*sphere*merged.data; 
[fullcomps_spec,fullcomps_freq] = spectopo(acts,0,merged.srate,'plot','off'); 
save('fullcomps_spec','fullcomps_spec'); save('fullcomps_freq','fullcomps_freq'); 
%}

fullcomps = load('fullcomps'); fullcomps = fullcomps.fullcomps; weights = fullcomps{1}; sphere = fullcomps{2}; 
acts = weights*sphere*merged.data; 

[pxx,f] = pwelch(acts',1000,250,480,merged.srate); 
allpxx(sb,:,:) = pxx; allwinv(sb,:,:) = pinv(weights*sphere); 

%{
fullcomps_spec = load('fullcomps_spec'); spec = fullcomps_spec.fullcomps_spec;
fullcomps_freq = load('fullcomps_freq'); freq = fullcomps_freq.fullcomps_freq; 
allspec(sb,:,:) = spec; 
allwinv(sb,:,:) = pinv(weights*sphere); 
%}
end


specindex = 120; 
normspec = (allpxx(:,:,:)); 
normspec = permute(normspec,[1,3,2]);
normspec = log(normspec); 

clear features; 

features(:,:,1:60) = mat2gray(normspec(:,:,1:60))-  repmat(mean(mat2gray(normspec(:,:,1:60)),2),[1,64,1]); 
%features(:,:,36:70) = mat2gray(diff(normspec(:,:,1:36),1,3)); 

%{
for i=1:64
   features(:,:,i+120) = normwinv(:,:,i); 
end
%}
res_features = reshape(features,[numel(features(:,:,1)),size(features,3)]); 
idx = kmeans(res_features,100,'Distance','correlation','Replicates',15); 
res_idx = reshape(idx,[8,64]); 
for i=1:100; disp(i); 
    figure,
    indsi = find(res_idx==i); 
    [x,y] = ind2sub(size(res_idx),indsi); 
    for j=1:length(indsi)
       subplot(ceil(sqrt(length(indsi))), ceil(sqrt(length(indsi))),j) ;
       topoplot(double(squeeze(allwinv(x(j),:,y(j)))),merged.chanlocs) ; 
    end
end

goods = [97,89,84,76,61,57,52,47,42,41,35,31,23,10,7,6,4,3]; 

good_comps = zeros(size(res_idx)); 
for i=1:length(goods)
    indsi = find(res_idx==goods(i)); 
    [x,y] = ind2sub(size(res_idx),indsi); 
    for j=1:length(x)
       good_comps(x,y) = 1;  
    end
end

for i=1:8
    goodinds = find(good_comps(i,:)==1); 
    figure,for j=1:length(goodinds) ; subplot(ceil(sqrt(length(goodinds))),ceil(sqrt(length(goodinds))),j) ; topoplot(squeeze(allwinv(i,:,goodinds(j))),merged.chanlocs); title(j); end
    cd(['E:\badger_eeg\',subs{i}]);
    fullcomps_goodinds = goodinds; save('fullcomps_goodinds','fullcomps_goodinds'); 
end


