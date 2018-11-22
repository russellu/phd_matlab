clear all ; close all; 

subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'};
goodcs = {[5,6,8,9,10,11,12,13,14,15,16,17,19,20,21,23,25,30,31,42,37],...
          [7,9,10,11,12,13,14,15,16,17,20,29,36,49,55],...
          [6,9,10,11,12,13,15,16,17,18,20,24,27,37,40,41],...
          [7,9,11,12,16,21,23,36],...
          [8,10,11,12,13,15,16,18,22,47],...
          [6,8,9,11,12,14,16,17,21,23,32,55],...
          [7,8,11,12,14,17,19,24,25,26,28,31,42,45,51,60],...
          [8,9,10,14,15,18,22,26,28,32,33,41]};

for sb=1:length(subs)
cd(['E:\badger_eeg\',subs{sb}]);
fullcomps = load('fullcomps'); fullcomps = fullcomps.fullcomps; weights = fullcomps{1}; sphere = fullcomps{2}; 
winv = pinv(weights*sphere); 
bads = zeros(1,64) ; bads(goodcs{sb}) = 1; bads = find(bads==0); 

denbcgs = dir('denbcg*gradeeg*set');
figure,
for denb=1:length(denbcgs)
    %if denb==1; merged = pop_loadset(denbcgs(denb).name); else merged = pop_mergeset(merged,pop_loadset(denbcgs(denb).name)); end
    eeg = pop_loadset(denbcgs(denb).name); 
    acts = weights*sphere*eeg.data; 
    acts(bads,:) = 0;
    invacts = (winv*acts); 
    subplot(2,3,denb) ; plot(eeg.data(46,:)) ; hold on ; plot(invacts(46,:)); 
    eeg.data = invacts; pop_saveset(eeg,['denica_',denbcgs(denb).name]);
end


%acts = weights*sphere*merged.data; 

%acts(bads,:) = 0;
%invacts = (winv*acts); 


%figure,plot(merged.data(46,:)) ; hold on ; plot(invacts(46,:)); title(subs{sb}); 


%{
[pxx,f] = (pwelch(acts',1000,250,480,merged.srate)); pxx = log(pxx); 

figure,
for i=1:64
   subplottight(10,14,i*2-1); topoplot(winv(:,i),merged.chanlocs); text(.55,0,num2str(i)); 
   tmax = max((pxx(1:60,i)));
   tmin = min((pxx(1:60,i)));    
   subplottight(10,14,i*2) ; plot(f(1:60),(pxx(1:60,i)),'LineWidth',2,'Color',[0,0,0]); vline(10,'r'); ylim([tmin,tmax]); set(gca,'XTick',[],'YTick',[]); 
   
end
%}

end





