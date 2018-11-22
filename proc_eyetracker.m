cd E:\russ_eye_date_07122018

rot2 = load_xdf('rotation_2.xdf'); 
ts1 = rot2{1};
ts2 = rot2{2}.time_series; 
ts3 = rot2{3}.time_series; 
ts4 = rot2{4}.time_series; 

pup1 = ts2(1,:); 
stamps1 = ts2(3,:);
pup2 = ts3(1,:); 
stamps2 = ts3(3,:); 

trigs = ts1.time_series; 
trigstamps = ts1.time_stamps; 
stimtrigs = {'S1','S2','S3','S4'};
for i=1:length(stimtrigs) ; lat_inds = find(strcmpi(stimtrigs{i},trigs)); stimlats(i,:) = trigstamps(lat_inds);  end

clear epochs; 
for i=1:size(stimlats,1)
    for j=1:size(stimlats,2)
        diff_ij = abs(stimlats(i,j) - stamps1); 
        ind_ij = find(diff_ij==min(diff_ij),1); 
        epochs(i,j,:) = pup1(ind_ij-120:ind_ij+60*20);       
    end
end


