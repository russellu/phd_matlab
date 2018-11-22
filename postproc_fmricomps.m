clear all ; close all 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','lisa','marc','marie','mathieu',...
    'maxime','mingham','patricia','po','russell','sunachakan','tah','thititip','vincent'} ; 

comps = {[18,13],[13,10],[10,8],[13,7],[12,23],[6,15],[18,8],[27,30],[16,20],[18,9],[10,9],[29,15],[6,8],[10,12],...
    [14,36],[13,8],[34,19],[8,16],[21,20],[9,13],[23,16],[17,12],[18,7],[11,25]} ;
for s=1:length(subs)
   cd(['C:\shared\all_white_normals\fmris\sub_',subs{s},'\melodic']) ; 
   compts = load('melodic_mix') ; 
   allpos(s,:) = compts(:,comps{s}(1)) ; 
   allneg(s,:) = compts(:,comps{s}(2)) ; 
end
allpos = eegfiltfft(allpos,0.5,0.01,1) ; 
allneg = eegfiltfft(allneg,0.5,0.01,1) ;

cd ../trigs ; 
stimTimes = load(['stimTimes',num2str(1),'.mat']) ; cd ..
stimTimes = stimTimes.stimTimes ; 
stimTimes = cell2mat(stimTimes) ; 
times = stimTimes(1:2:end) ; trs = round(times/2) ; 

meanpos = squeeze(mean(allpos,1)) ; 
meanneg = squeeze(mean(allneg,1)) ; 
save('meanpos','meanpos') ; 
save('meanneg','meanneg') ; 


shadedErrorBar([],mean(allpos,1),std(allpos,0,1)/sqrt(24),'r') ; 
hold on ; shadedErrorBar([],mean(allneg,1),std(allneg,0,1)/sqrt(24),{'b'}) ;
vline(trs,'k') ; 

