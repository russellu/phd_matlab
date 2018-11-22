cd E:\data_for_Russ
names = {'eyes_closed_01.vhdr','eyes_open_01.vhdr'};
savenames = {'ec_01','eo_01'};
ts = {32582:32592,58446:58456};
for st=1:length(names)
eeg  = pop_loadbv('.',names{st});

rawts = eeg.data; 
for i=1:64
    rawt_i = eeg.data(i,:); 
    template_i = eeg.data(i,ts{st});
    conved_i = conv(rawt_i,template_i,'same');
    rawts(i,:) = conved_i(1:size(eeg.data,2)); 
end
mrawts = mean(rawts,1)./max(mean(rawts,1)); 
mrawts(mrawts<0.3) = 0; 
[pks,locs] = findpeaks(double(mrawts),'MinPeakDistance',eeg.srate*3); 
figure,plot(eeg.data(3,:)) ; hold on ; vline(locs)
types = {}; newlocs = {}; 
for i=1:length(locs); types{i} = 'TMS'; newlocs{i} = locs(i);  end; 
eeg = eeg_addnewevents(eeg,newlocs,types); 

eps = pop_epoch(eeg,{'TMS'},[-3,3]); 
%plot(squeeze(eps.data(3,:,:)))
pop_saveset(eeg,savenames{st}); 

end




