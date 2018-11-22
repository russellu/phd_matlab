clear all ; close all; 

subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'};
hrfacs = [3.1,3.5,3.5,3.5,3.5,4,3.5,3.5];

for sb=1:length(subs)
cd(['E:\badger_eeg\',subs{sb}])
bades = [61,29,62,30,55,11,53,1,63,2,54,12,56,14,58,16,60,10,20,9,59,15,57,13]; 

sets = {'gradeeg*allstim*01*set','gradeeg*allstim*02*set','gradeeg*gamma*01*set','gradeeg*gamma*02*set','gradeeg*movie*set','gradeeg*rest*set'}; 
for st=1:length(sets)
    setname = dir(sets{st}); 
    eeg = pop_loadset(setname.name); 
    hrchan = eeg.data(32,:); 
    eeg.data(32,:) = rand(1,length(eeg.data))/10; 
    filtdat = eegfiltfft(eeg.data,eeg.srate,1,128); 
    [weights,sphere] = runica(filtdat(:,1:4:end),'maxsteps',128); 
    winv = pinv(weights*sphere); 
    acts = weights*sphere*filtdat; 
    clear skews; 
    for i=1:10
        icount =1 ; 
        for j=1:500:length(acts)-1000
            skews(i,icount) = skewness(acts(i,j:j+300)); icount = icount + 1; 
        end
    end
    
    [sv,si] = sort(abs(median(skews,2)),'descend'); 
     %figure,for i=1:4 ; subplot(4,1,i); plot(acts(si(i),5000:15000)*median(skews(si(i),:),2)); end
    
     xc1 = xcorr(acts(si(1),:),300,'coeff');
     xc2 = xcorr(acts(si(2),:),300,'coeff'); 
     
     mxc = xc1+xc2; mxc = mxc(300:end) ; mxc(1:150) = 0; 
     hrate = find(mxc==max(mxc)); 
     
     [pks1,locs1] = findpeaks(acts(si(1),:)*median(skews(si(1),:),2),'MinPeakDistance',hrate-hrate/hrfacs(sb)); 
     [pks2,locs2] = findpeaks(acts(si(2),:)*median(skews(si(2),:),2),'MinPeakDistance',hrate-hrate/hrfacs(sb)); 

     
     figure, plot(acts(si(1),:)*median(skews(si(1),:),2)) ; vline(locs1); title(subs{sb}); 
     
     neweeg = eeg;     si(11:64) = 11:64; 
     neweeg.data = acts(si,:);
     if st==1 ; merged = neweeg; else merged = pop_mergeset(merged,neweeg); end
     
     acts(6:end,:) = 0; 
     bcgeeg = eeg;
     bcgeeg.data = winv*acts; 
     if st==1 ; bcgmerged = bcgeeg; else bcgmerged = pop_mergeset(bcgmerged,bcgeeg); end

    part1 = floor(hrate*.4); part2 = floor(hrate*0.7); 

    clear peakinds peakvals; 
    icount=1; 
    for i=2:length(locs1)-1
        peakinds(icount,:) = locs1(i)-part1:locs1(i)+part2; 
        peakvals(icount,:,:) = filtdat(:,locs1(i)-part1:locs1(i)+part2); 
        icount = icount + 1; 
    end

    %figure,imagesc(squeeze(peakvals(:,46,:)),[-100,100]); colormap jet; 
    
    clear bdiffs; 
    for i=1:length(bades)
       badchani = squeeze(peakvals(:,bades(i),:));  
       for j=1:size(badchani,1)
          diffs_ij = repmat(badchani(j,:),[size(badchani,1),1]); 
          diffmat = (diffs_ij - badchani);       
          bdiffs(i,j,:) = sum(abs(diffmat),2); 

       end
    end

    mbdiffs = squeeze(mean(bdiffs,1)); 
    [sv,si] = sort(mbdiffs,2,'ascend') ;
    templates = zeros(size(peakvals)); 
    for i=1:size(si,1)
        templates(i,:,:) = squeeze(mean(peakvals(si(i,2:25),:,:),1));     
    end

    neweeg = eeg ; 
    for i=1:size(peakinds,1)
       neweeg.data(:,peakinds(i,:)) = eeg.data(:,peakinds(i,:)) - squeeze(templates(i,:,:)); % squeeze(subbed(i,:,:));      
    end

    figure,plot(eeg.data(46,:)) ; hold on ; plot(neweeg.data(46,:)); 
    pop_saveset(neweeg,['denbcg_',strrep(setname.name,'.vhdr','.set')]);

    if st==1; denbcg_merged = neweeg; else denbcg_merged = pop_mergeset(denbcg_merged,neweeg); end
    
    
    %{
    hrxcs = xcorr(hrchan,350,'coeff'); 
    for i=1:10 ; xcs(i,:) = xcorr(acts(i,:),350,'coeff'); end
    
    corrs = corr(hrxcs(80:300)',(xcs(:,80:300))') ;
    [sv,si] = sort(corrs,'descend'); 
    figure,for i=1:4 ; subplot(4,1,i); plot(acts(si(i),5000:15000)); end
    %}
    
    %{
    eeg2 = pop_loadbv('.',setname.name); 
    reseeg2 = pop_resample(eeg2,250); 
    gradeeg = remove_gradient2(eeg2); 
    gradeeg = pop_chanedit(gradeeg,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
    pop_saveset(gradeeg,['gradeeg_',strrep(setname.name,'.vhdr','.set')]);
    %}
    
end
pop_saveset(denbcg_merged,'denbcg_merged.set'); 
%pop_saveset(bcgmerged,'bcgmerged.set');
%pop_saveset(merged,'mergedacts.set'); 
end




%{
[weights,sphere] = runica(neweeg,'maxsteps',128); 
winv = pinv(weights*sphere); 
figure,for i=1:64 ; subplot(5,13,i); topoplot(winv(:,i),gradeeg.chanlocs) ; end
acts = weights*sphere*neweeg; 

[spec,hz] = spectopo(acts,0,gradeeg.srate,'plot','off'); 
newgrad = gradeeg; newgrad.data = acts; 
allep = pop_epoch(newgrad,{'S  1','S  2'},[-1,6]); clear mstersp; 
 for j=1:64  ; disp(j); 
            [mstersp(j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(j,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',80,'baseline',0,'verbose','off','timesout',200) ; 
 end

figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mstersp(i,:,:)),[-6,6]) ; axis xy ; colormap jet; end
%}
 
 

