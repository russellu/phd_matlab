clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','valerie','tegan'} ; 
% component labels
categories = {'bcg','center','eyes','rim'... %
              'motor_left','motor_center','motor_right',... % 
              'frontal_left','frontal_center','frontal_right',... % 7
              'occ_left','occ_center','occ_right',... % 
              'parietal_left','parietal_center','parietal_right',... % 
              'temporal_left','temporal_right'... %
              'noise_symmetric','noise_focal','noise_broad'} ; % 
          
alex = {'bcg','bcg','bcg','center','motor_right','noise_broad','frontal_center','frontal_left','eyes','temporal_left','occ_left','parietal_left','occ_right','parietal_center','occ_right','frontal_right','occ_left',...
        'parietal_center','parietal_left','temporal_left','temporal_right','occ_center','eyes','noise_broad','noise_broad','occ_right','parietal_right','noise_broad','noise_symmetric','noise_broad','noise_broad',...
        'noise_broad','noise_broad','noise_focal','noise_broad','noise_broad','noise_focal','noise_broad','noise_broad','noise_focal','eyes','noise_focal','eyes','noise_focal','noise_broad','noise_focal','noise_broad',...
        'noise_broad','noise_broad','noise_focal','noise_focal','noise_focal','noise_broad','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal',...
        'noise_focal','noise_focal'}; 
    
dina = {'bcg','bcg','bcg','bcg','center','occ_right','bcg','motor_right','frontal_center','eyes','rim','motor_center','motor_left','occ_left','eyes','occ_right','parietal_left','frontal_right','bcg',...
        'temporal_left','noise_broad','noise_broad','noise_focal','noise_broad','noise_broad','noise_broad','occ_right','noise_broad','noise_broad','eyes','noise_broad','temporal_right','noise_focal','noise_focal',...
        'noise_broad','noise_broad','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_broad','noise_focal','noise_focal','noise_focal','noise_broad','noise_focal','noise_focal','noise_broad',...
        'noise_broad','noise_focal','noise_focal','bcg','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal'} ; 
    
genevieve = {'bcg','bcg','bcg','center','eyes','frontal_right','occ_left','temporal_right','frontal_left','parietal_center','occ_right','motor_center','occ_center','temporal_left','noise_broad','occ_left','eyes',...
            'noise_broad','occ_right','temporal_right','noise_broad','motor_right','parietal_left','frontal_left','frontal_left','eyes','eyes','eyes','noise_focal','occ_right','noise_focal','noise_broad',...
            'noise_focal','noise_focal','noise_symmetric','temporal_right','noise_broad','noise_focal','parietal_left','noise_broad','noise_broad','noise_focal','noise_broad','noise_broad','eyes','noise_broad','noise_focal',...
            'noise_broad','noise_broad','noise_focal','noise_broad','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_broad','noise_focal','noise_focal','noise_focal',...
            'noise_focal','noise_focal'} ; 
        
jeremie = {'bcg','bcg','bcg','center','eyes','bcg','occ_center','frontal_right','motor_right','bcg','parietal_center','noise_symmetric','noise_broad','bcg','eyes','noise_broad','temporal_left','noise_symmetric',...
    'noise_symmetric','parietal_left','noise_broad','occ_center','noise_focal','frontal_right','noise_broad','motor_left','frontal_right','noise_focal','noise_focal','parietal_right','temporal_right','noise_focal','noise_broad'...
    'noise_focal','noise_broad','noise_focal','noise_focal','noise_broad','noise_broad','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','eyes','eyes','noise_focal','noise_focal','noise_broad','noise_focal'...
    'noise_focal','noise_broad','noise_focal','noise_broad','noise_broad','noise_focal','noise_focal','noise_focal','noise_focal','noise_broad','noise_focal','noise_focal','noise_focal','noise_focal'} ; 

karl = {'bcg','bcg','bcg','center','bcg','bcg','frontal_center','bcg','occ_center','parietal_center','noise_Broad','eyes','bcg','bcg','noise_broad','bcg','eyes','occ_right','temporal_left','parietal_left','noise_broad'...
    'motor_right','noise_symmetric','occ_center','temporal_right','noise_broad','noise_broad','bcg','noise_broad','noise_broad','noise_broad','noise_broad','noise_broad','noise_broad','noise_broad','noise_focal','noise_broad',...
    'noise_broad','noise_broad','noise_broad','eyes','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','frontal_right','noise_focal','noise_broad','noise_focal','noise_focal','noise_broad','noise_broad',...
    'noise_focal','noise_Broad','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_broad','noise_focal','noise_broad','noise_focal'} ; 

russell = {'bcg','bcg','bcg','bcg','bcg','eyes','bcg','center','noise_broad','rim','motor_right','eyes','parietal_center','motor_center','noise_broad','parietal_left','occ_right','noise_broad','frontal_right','temporal_left',...
    'noise_broad','noise_broad','occ_left','noise_broad','frontal_left','noise_broad','noise_broad','noise_broad','noise_broad','noise_broad','noise_broad','noise_broad','noise_broad','eyes','noise_broad','noise_broad',...
    'noise_broad','noise_broad','noise_broad','noise_focal','noise_broad','noise_focal','noise_broad','noise_focal','noise_focal','noise_focal','noise_focal','noise_broad','noise_focal','noise_broad','noise_broad',...
    'noise_broad','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_broad','noise_broad','noise_focal','noise_focal','noise_focal','noise_focal'} ; 

sukhman = {'bcg','center','bcg','eyes','frontal_center','motor_center','occ_center','bcg','bcg','parietal_left','parietal_center','noise_broad','eyes','parietal_right','noise_broad','occ_center','frontal_right','noise_symmetric',...
           'noise_broad','noise_broad','parietal_center','temporal_left','noise_broad','occ_center','noise_broad','noise_broad','noise_broad','noise_broad','noise_broad','motor_right','noise_broad','noise_broad','noise_broad',...
           'noise_broad','noise_broad','noise_broad','noise_broad','noise_focal','noise_broad','noise_focal','noise_broad','noise_broad','noise_broad','noise_broad','noise_broad','noise_broad','noise_broad','noise_focal',...
           'noise_broad','noise_focal','noise_focal','noise_focal','noise_broad','noise_broad','noise_focal','noise_focal','noise_focal','noise_broad','noise_broad','noise_Broad','noise_focal','noise_broad','noise_broad',...
           'noise_broad'} ; 
       
valerie = {'bcg','bcg','bcg','bcg','center','eyes','occ_center','bcg','noise_broad','noise_broad','noise_broad','frontal_center','motor_center','rim','temporal_left','noise_broad','occ_center','noise_broad','noise_broad',...
    'noise_broad','noise_broad','parietal_center','noise_broad','eyes','motor_center','noise_symmetric','noise_broad','bcg','bcg','noise_focal','noise_broad','noise_broad','parietal_left','noise_focal','noise_broad',...
    'noise_focal','occ_center','noise_broad','noise_focal','frontal_left','noise_focal','noise_broad','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal',...
    'noise_focal','noise_focal','eyes','noise_focal','noise_broad','noise_focal','noise_broad','noise_focal','noise_broad','noise_focal','noise_focal','noise_focal','noise_focal'} ; 

tegan = {'bcg','bcg','eyes','center','bcg','bcg','frontal_center','occ_left','eyes','parietal_left','noise_broad','noise_broad','noise_focal','noise_symmetric','noise_broad','occ_right','parietal_right','frontal_center',...
    'noise_broad','noise_focal','frontal_left','frontal_left','occ_center','occ_right','noise_broad','noise_broad','noise_broad','parietal_center','parietal_center','temporal_left','noise_broad','noise_focal','noise_broad',...
    'noise_broad','noise_broad','occ_left','noise_broad','noise_focal','eyes','noise_broad','noise_Broad','noise_broad','noise_focal','noise_focal','noise_broad','noise_broad','noise_broad','noise_broad','noise_focal',...
    'noise_focal','noise_focal','noise_focal','noise_broad','noise_focal','noise_focal','noise_broad','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal','noise_focal'} ;

complabs = {alex,dina,genevieve,jeremie,karl,russell,sukhman,valerie,tegan} ;

for sub=1:length(subs)
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    prefix = 'allfreq_' ; 
    rests=dir('allfreq_*rest*set') ;     
    allrest{sub} = pop_loadset(rests(1).name) ; 
end
EEG = allrest{1} ; 

% plot the topoplots
elocs = {allrest{1}.chanlocs.labels} ; 
for s=1:9
    figure,
    for i=1:64 ; subplot(5,13,i) ; 
        topoplot(allrest{1}.icawinv(:,i),allrest{1}.chanlocs) ; title(i) ; 
    end    
    allkmaps(s,:,:) = allrest{s}.icawinv(:,:) ; 
    suptitle(subs{s}) ; 
end

% test the individual subjects
for i=1:64 ; subplot(5,13,i) ; 
    topoplot(allrest{9}.icawinv(:,i),allrest{1}.chanlocs) ; title(tegan{i}) ;  
end    

% match the components to the list, and sort them
catcount = ones(1,length(categories)) ; 
alltopos = cell(1,length(categories)) ; 
for i=1:length(complabs) 
    subi = complabs{i} ; 
    for j=1:length(subi)
        labij = subi{j} ; 
        catind = find(strcmpi(labij,categories)) ;
        vec = allrest{i}.icawinv(:,j) ; 
        alltopos{catind}{catcount(catind)} = vec ; 
        allcompinds(i,j) = catind ; 
        catcount(catind) = catcount(catind) + 1 ; 
    end
end

ind = 12 ; figure,
for i=1:length(alltopos{ind}) ;
   subplot(5,6,i) ;
   topoplot(squeeze(alltopos{ind}{i}),EEG.chanlocs) ; 
    
end

figure,
for i=1:length(alltopos)
   subplot(3,7,i) ;  clear topos
   for j=1:length(alltopos{i})
      topos(j,:) = squeeze(alltopos{i}{j}) ;  
   end
   topoplot(squeeze(mean((topos),1)),EEG.chanlocs) ; title(i) ; 

end


for s=1:length(subs)
    cd(['c:/shared/badger_eeg/',subs{s}]) ; 
    inds_s = allcompinds(s,:) ; 
    save('compinds','inds_s') ; 
    
    
    
end









