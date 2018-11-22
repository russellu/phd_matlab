clear all ; close all; 
cd E:\russ_eye_date_07122018\csv_data;
triggers = {'S2','S4','S12','S14','S22','S24','S32','S34','S42','S44'};
eye = 'eye0';

for i=1:length(triggers)
   csv_i = dir(['*',eye,'_data_original_trigger_',triggers{i},'.csv']); 
    
   data = load(csv_i.name); 
  size(data)
   % alldata(i,:,:) = data; 
    
end

