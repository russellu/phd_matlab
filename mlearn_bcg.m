clear all ; close all; 

subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','tegan','valerie'};
for sb=9%:length(subs)
cd(['E:\badger_eeg\',subs{sb}])

sets = {'*allstim*01*vhdr','*allstim*02*vhdr','*gamma*01*vhdr','*gamma*02*vhdr','*movie*vhdr','*rest*vhdr'}; 
for st=1:length(sets)
    setname = dir(sets{st}); 
eeg2 = pop_loadbv('.',setname.name); 
reseeg2 = pop_resample(eeg2,250); 

gradeeg = remove_gradient2(eeg2); 
gradeeg = pop_chanedit(gradeeg,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;

filt = eegfiltfft(gradeeg.data,gradeeg.srate,1,128); % high pass
bades = [61,29,62,30,55,11,53,1,63,2,54,12,56,14,58,16,60,10,20,9,59,15,57,13]; 

crosscorrs = xcorr(smooth(double(gradeeg.data(32,:)),25),500); 
xcorrs = crosscorrs(500:end); xcorrs(1:170) = 0; maxhr = find(xcorrs==max(xcorrs)); 


[pks,locs] = findpeaks(smooth(double(gradeeg.data(32,:))),'MinPeakDistance',maxhr-maxhr/6); 

part1 = floor(maxhr*.2); part2 = floor(maxhr*0.9); 

clear peakinds peakvals; 
icount=1; 
for i=2:length(locs)-1
    peakinds(icount,:) = locs(i)-part1:locs(i)+part2; 
    peakvals(icount,:,:) = filt(:,locs(i)-part1:locs(i)+part2); 
    icount = icount + 1; 
end

figure,imagesc(squeeze(peakvals(:,46,:)),[-100,100]); colormap jet; 

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

%subbed = peakvals - templates; 

neweeg = gradeeg ; 
for i=1:size(peakinds,1)
   neweeg.data(:,peakinds(i,:)) = gradeeg.data(:,peakinds(i,:)) - squeeze(templates(i,:,:)); % squeeze(subbed(i,:,:));      
end

figure,plot(gradeeg.data(46,:)) ; hold on ; plot(neweeg.data(46,:)); 

pop_saveset(neweeg,['denbcg_',strrep(setname.name,'.vhdr','.set')]);

end

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
 
 

