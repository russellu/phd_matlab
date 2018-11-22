clear all ; close all;
subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'};

comps = {[5,23,21],[8,21,20],[5,6,15],[5,13,19],[5,9,14],[4,11,18],[3,9,13],[8,13,16],[7,9,14]}; 

amps = zeros(9,1,3,64,50,15); 
stds = zeros(9,1,3,64,50,15); 
for sb=1:length(subs)
    disp(sb); 
    cd(['E:\badger_eeg\',subs{sb},'\outside']);
    
    stims = {'S  2'};
    basetime = [-500,0];
    tasktime = [500,3000]; 
    posttime = [3500,4000]; 

    hzs = 2:2:100; 
    for hz=1:length(hzs) ; disp(subs{sb}); 
       eeg = pop_loadset(['gen_hz_',num2str(hzs(hz)),'.set']); 
       for st=1:length(stims)
           eps = pop_epoch(eeg,{stims{st}},[-1,4]); 
           eps.data = abs(eps.data); 
           amps(sb,st,1,:,hz,:) = squeeze(mean(eps.data(:,(eps.times > basetime(1) & eps.times< basetime(2)),1:end),2)); 
           amps(sb,st,2,:,hz,:) = squeeze(mean(eps.data(:,(eps.times > tasktime(1) & eps.times< tasktime(2)),1:end),2)); 
           amps(sb,st,3,:,hz,:) = squeeze(mean(eps.data(:,(eps.times > posttime(1) & eps.times< posttime(2)),1:end),2)); 
           stds(sb,st,1,:,hz,:) = squeeze(std(eps.data(:,(eps.times > basetime(1) & eps.times< basetime(2)),1:end),0,2)); 
           stds(sb,st,2,:,hz,:) = squeeze(std(eps.data(:,(eps.times > tasktime(1) & eps.times< tasktime(2)),1:end),0,2)); 
           stds(sb,st,3,:,hz,:) = squeeze(std(eps.data(:,(eps.times > posttime(1) & eps.times< posttime(2)),1:end),0,2)); 
       end
    end
end
% clean the bad trials:

ampdiff = squeeze(amps(:,:,2,:,:,:) - amps(:,:,1,:,:,:)) ;
allampdiff = ampdiff; save('allampdiff','allampdiff'); 
allstds = stds; save('allstds','allstds'); 
substds = (squeeze(mean(mean(mean(stds(:,:,1:2,:,20:end,:),3),4),5)));
[sv,si] = sort(substds,3,'descend'); 

mampdiff = zeros(size(ampdiff,1),size(ampdiff,2),size(ampdiff,3)); 
for i=1:9
        mampdiff(i,:,:) = squeeze(mean(ampdiff(i,:,:,:),4)); 

end
mmampdiff = squeeze(mean(mampdiff,2)); 
postchans = [60,61,62,63,64,29,30,31,23,56,24,57,25,58,26,59,27];
c = corr(squeeze(mean(mmampdiff(:,postchans,:),2)));
save('mampdiff','mampdiff');














%{
for sb=1:length(subs)
    disp(sb); 
    cd(['E:\badger_eeg\',subs{sb},'\outside']);
    efile = dir('*vhdr');
    eeg = pop_loadbv('.',efile(1).name); 
    eeg = pop_resample(eeg,250); 
    eeg.data = eeg.data - eegfiltfft(eeg.data,eeg.srate,59,61); 
    eeg = pop_chanedit(eeg,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
    cleanweights = load('compdat'); cleanweights = cleanweights.compdat; 
    weights = cleanweights{1}; sphere = cleanweights{2}; 
    winv = pinv(weights*sphere); 
    neweeg = eeg;
    acts = weights*sphere*eeg.data;     
    invacts = winv(:,comps{sb})*acts(comps{sb},:); 
    neweeg.data = invacts; 
    
    freqs = 1:2:120;
    clear freqtopos; 
    %invacts = invacts - eegfiltfft(invacts,eeg.srate,59,61) - eegfiltfft(invacts,eeg.srate,84,86); 
    for f=1:length(freqs)    
        gammafilt = eegfiltfft(invacts,eeg.srate,freqs(f)-1,freqs(f)+1) ;
        gammaeeg = eeg; gammaeeg.data = gammafilt; 
        pop_saveset(gammaeeg,['gen_hz_',num2str(f*2),'.set']);    
    end
    
    gammafilt = eegfiltfft(invacts,eeg.srate,8,25) ;
    gammaeeg = eeg; gammaeeg.data = gammafilt; 
    pop_saveset(gammaeeg,['gen_alpha_hz.set']);   
    gammafilt = eegfiltfft(invacts,eeg.srate,40,90) ;
    gammaeeg = eeg; gammaeeg.data = gammafilt; 
    pop_saveset(gammaeeg,['gen_gamma_hz.set']);   
       
end

%}





%{
amps = zeros(9,1,3,64,50,130); 
stds = zeros(9,1,3,64,50,130); 
for sb=1:length(subs)
    disp(sb); 
    cd(['E:\badger_eeg\',subs{sb},'\outside']);
    
    stims = {'S  2'};
    basetime = [-500,0];
    tasktime = [500,3000]; 
    posttime = [3500,4000]; 

    hzs = 2:2:100; 
    for hz=1:length(hzs) ; disp(subs{sb}); 
       eeg = pop_loadset(['gen_hz_',num2str(hzs(hz)),'.set']); 
       for st=1:length(stims)
           eps = pop_epoch(eeg,{stims{st}},[-1,4]); 
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
%}

%{
for sb=1:length(subs)
    
    cd(['E:\badger_eeg\',subs{sb},'\outside']);
    ls
    efile = dir('*vhdr');
    eeg = pop_loadbv('.',efile(1).name); 
    eeg = pop_resample(eeg,250); 
    eeg.data = eeg.data - eegfiltfft(eeg.data,eeg.srate,59,61); 
    eeg = pop_chanedit(eeg,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;

    filtdat = eegfiltfft(eeg.data,eeg.srate,1,128);
    [weights,sphere] = runica(filtdat,'maxsteps',128); 
    compdat{1} = weights; compdat{2} = sphere; 
    save('compdat','compdat'); 
    winv = pinv(weights*sphere); 
    neweeg = eeg; neweeg.data = weights*sphere*eeg.data; 
    
    trigs = {'S  2'}; 
    clear ersp ; 
    epica = pop_epoch(neweeg,{trigs{1}},[-1,6]); 
     for i=1:64
            [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(i,:,:)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',0,'verbose','off','timesout',100) ; 
     end
    
    figure,for i=1:64 ; subplot(8,16,i*2) ; imagesc(squeeze(ersp(i,:,:)),[-5,5]) ; title((i)); axis xy ; colormap jet; end
    for i=1:64 ; subplot(8,16,i*2-1) ; topoplot(winv(:,(i)),eeg.chanlocs); end
    
    
end
%}