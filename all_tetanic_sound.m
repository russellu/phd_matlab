clear all ; close all ; 
stims = {'S  1','S  2','S 11','S 12'} ;
subs = {'AD','AM','CV','DT','FMP','KC','LB','RB','TD'};
bades = {[],[],[21,61],[9,42],[19,57,62],[1,32,33,37],[10],[],[18,50]};

cd(['E:\jf_data\Tetanic Visual\RB']) ; rmerged = pop_loadset('merged.set'); 

components = {25:35,40:48,56:64,70:76,85:91};

for sb=1:length(subs)
    
    cd(['E:\jf_data\Tetanic Sound\',subs{sb}]) ;
    %{
    
    discr = dir('*vhdr') ;
    for i=1:length(discr)
       EEG = pop_loadbv('.',discr(i).name) ; 
       EEG = pop_resample(EEG,256) ; 
       if i==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG) ; end
    end
    merged = pop_chanedit(merged,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
    merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61) - eegfiltfft(merged.data,merged.srate,84,86) - eegfiltfft(merged.data,merged.srate,0,1);  

    %figure,bar(sum(abs(diff(merged.data,1,2)),2))
    %merged = pop_interp(merged,[32],'spherical');
    pop_saveset(merged,'merged.set'); 
    %}
    
    merged = pop_loadset('merged.set');   
   
    for i=1:length(bades{sb})
       merged.data(bades{sb}(i),:) = rand(1,length(merged.data))/100;  
        
    end
    
    merged.data = eegfiltfft(merged.data,merged.srate,1,128); 
    
    eps = pop_epoch(merged,stims,[-1,2]);
    %figure; for i=1:64 ; subplot(5,13,i) ; plot(zscore(squeeze(mean(eps.data(i,300:450,:),3)))) ; ylim([-7,7]); title(i); end ; suptitle(subs{sb}); 
    
    
    
    for i=1:64
    alleps(sb,i,:) = (squeeze(mean(eps.data(i,300:450,:),3))); 
    end
    
    elecdiffs = squeeze(mean(mean(eps.data(:,[components{3},components{5}],:),2),3)-mean(mean(eps.data(:,[components{2},components{4}],:),2),3)); 
    [sv,si] = sort(elecdiffs,'descend'); 
    
    elec_inds = si ; save('elec_inds','elec_inds'); 
    
    for i=1:length(stims)
       eps_i = pop_epoch(merged,{stims{i}},[-1,2]);  
       
       zeps = (eps_i.data); 
        mdiffs(sb,i,:) = -squeeze(mean(mean(zeps(elec_inds(1:10),:,:),3),1)); 
    end
    
    
    
    
    
    
end

colors = {'g','c','m','y','k'};


times = 300:450; 
etimes = eps_i.times(times); 

subplot(3,5,1); 
plot(eps_i.times(times),squeeze(mean(mdiffs(:,1,times),1)),'r','LineWidth',2);  hold on ;
plot(eps_i.times(times),squeeze(mean(mdiffs(:,3,times),1)),'b','LineWidth',2); ylabel('\muV'); xlabel('time(ms)'); 
for i=1:5 ; vline(etimes(components{i}),colors{i}) ; end 
 
subplot(3,5,2); 
plot(eps_i.times(times),squeeze(mean(mdiffs(:,2,times),1)),'r','LineWidth',2);  hold on ;
plot(eps_i.times(times),squeeze(mean(mdiffs(:,4,times),1)),'b','LineWidth',2); ylabel('\muV'); xlabel('time(ms)');
for i=1:5 ; vline(etimes(components{i}),colors{i}) ; end 
titles = {'P1','N1','P2','N2','P3'}; 

subplot(3,5,3) ; for i=1:5 ; plot(i,colors{i}) ; hold on  ;end ; legend(titles); 

tmdiffs = mdiffs(:,:,times); 
for i=1:5
   subplot(3,10,i+10); 
    barwitherr(squeeze(std(mean(tmdiffs(:,[1,3],components{i}),3),0,1))/3,squeeze(mean(mean(tmdiffs(:,[1,3],components{i}),3),1)));
    [h,p,ci,stats] = ttest(squeeze((mean(tmdiffs(:,[3],components{i}),3))),squeeze((mean(tmdiffs(:,[1],components{i}),3))));
    title([titles{i},': ',format_t(stats.tstat),' ',format_p(p)]); set(gca,'XTickLabel',{'pre','post'}); ylabel('\muV'); 
end

for i=1:5
   subplot(3,10,i+20); 
    barwitherr(squeeze(std(mean(tmdiffs(:,[2,4],components{i}),3),0,1))/3,squeeze(mean(mean(tmdiffs(:,[2,4],components{i}),3),1)));
    [h,p,ci,stats] = ttest(squeeze((mean(tmdiffs(:,[4],components{i}),3))),squeeze((mean(tmdiffs(:,[2],components{i}),3))));
    title([titles{i},': ',format_t(stats.tstat),' ',format_p(p)]); set(gca,'XTickLabel',{'pre','post'});ylabel('\muV'); 
end



