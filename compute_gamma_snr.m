clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 
amps = zeros(24,6,3,64,50,130); 
stds = zeros(24,6,3,64,50,130); 
for sb=1:length(subs)
    disp(sb); 
    cd(['E:\clean_allres\',subs{sb}]) ; ls 
    
    stims = {'S 11','S 12','S 13','S 14','S 15','S 16'};
    basetime = [-500,0];
    tasktime = [500,2000]; 
    posttime = [2500,3000]; 

    hzs = 2:2:100; 
    for hz=1:length(hzs) ; disp(subs{sb}); 
       eeg = pop_loadset(['gen_hz_',num2str(hzs(hz)),'.set']); 
       for st=1:length(stims)
           eps = pop_epoch(eeg,{stims{st}},[-1,3]); 
           eps.data = abs(eps.data); 
           amps(sb,st,1,:,hz,:) = squeeze(mean(eps.data(:,(eps.times > basetime(1) & eps.times< basetime(2)),1:130),2)); 
           amps(sb,st,2,:,hz,:) = squeeze(mean(eps.data(:,(eps.times > tasktime(1) & eps.times< tasktime(2)),1:130),2)); 
           amps(sb,st,3,:,hz,:) = squeeze(mean(eps.data(:,(eps.times > posttime(1) & eps.times< posttime(2)),1:130),2)); 
           stds(sb,st,1,:,hz,:) = squeeze(std(eps.data(:,(eps.times > basetime(1) & eps.times< basetime(2)),1:130),0,2)); 
           stds(sb,st,2,:,hz,:) = squeeze(std(eps.data(:,(eps.times > tasktime(1) & eps.times< tasktime(2)),1:130),0,2)); 
           stds(sb,st,3,:,hz,:) = squeeze(std(eps.data(:,(eps.times > posttime(1) & eps.times< posttime(2)),1:130),0,2)); 
       end
    end
end
% clean the bad trials:

ampdiff = squeeze(amps(:,:,2,:,:,:) - amps(:,:,1,:,:,:)) ;
allampdiff = ampdiff; save('allampdiff','allampdiff'); 
allstds = stds; save('allstds','allstds'); 
substds = (squeeze(mean(mean(mean(stds(:,:,1:2,:,20:end,:),3),4),5)));
[sv,si] = sort(substds,3,'descend'); 

mampdiff = zeros(size(ampdiff,1),size(ampdiff,2),size(ampdiff,3),size(ampdiff,4)); 
for i=1:24
    for j=1:6
        mampdiff(i,j,:,:) = squeeze(mean(ampdiff(i,j,:,:,si(i,10:end)),5)); 
    end
end
mmampdiff = squeeze(mean(mampdiff,2)); 
postchans = [60,61,62,63,64,29,30,31,23,56,24,57,25,58,26,59,27];
c = corr(squeeze(mean(mmampdiff(:,postchans,:),2)));
save('mampdiff','mampdiff');





