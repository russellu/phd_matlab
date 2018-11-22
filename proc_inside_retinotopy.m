
clear all  ; close all ;


comps = {[31,19],[17,18],[7,30],[8,16],[27,20],1:64} ; % left gamma first
subs = {'alex','dina','genevieve','karl','russell','tegan'} ; 

for sub=6%:length(subs) ; 
cd(['c:/shared/badger_eeg/',subs{sub}]) ; ls  ;
sounds=dir('highfreq*allstim*set') ;
for i=1:max(size(sounds)) ; 
   EEG = pop_loadset(sounds(i).name) ; 
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
   alleegs{i} = EEG ; 
end
trigs{1} = {'S 11','S 12','S 13','S 14'} ; trigs{2} = {'S 21','S 22','S 23','S 24'} ; trigs{3} = {'S 31','S 32','S 33','S 34'} ; trigs{4} = {'S 41','S 42','S 43','S 44'} ; 
trigs{5} = {'S 51','S 52','S 53','S 54'} ; trigs{6} = {'S 61','S 62','S 63','S 64'} ; trigs{7} = {'S 71','S 72','S 73','S 74'} ; trigs{8} = {'S 81','S 82','S 83','S 84'} ; 
clear ersp ; 
for i=1:2
    for t=1:length(trigs)
    ep = pop_epoch(alleegs{i},trigs{t},[-1,12]) ; 
        for c=1:64 ; 
            [ersp(i,t,c,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(c,:,:)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                        'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'timesout',150) ; 
        end
    end
end
mersp = squeeze(mean(ersp,1)) ; 
startangles = mod(210 - [0,60,120,180,240,300],360) ;
tpoints = find(times>0 & times<10) ; 
angincr = 360/length(tpoints) ; 
clear angleinds
for i=1:length(startangles) ; angleinds(i,:) = round(mod(startangles(i):-angincr:startangles(i)-360+1,360)) ; end
uniques = unique(angleinds) ; offset = length(find(times<0)) ; 
anglepower = zeros(size(mersp,2),size(mersp,3),size(angleinds,2)) ; 
for i=1:size(angleinds,1) 
    anglepoweri = squeeze(mersp(i,:,:,tpoints)) ; 
    [sv,si] = sort(angleinds(i,:),'descend') ; % sort the wedge angles
    anglepower = anglepower + anglepoweri(:,:,si) ; 
end
anglepower = anglepower/6 ; 

lefts(sub,:,:) = squeeze(anglepower(comps{sub}(1),:,:)) ; 
rights(sub,:,:) = squeeze(anglepower(comps{sub}(2),:,:)) ; 
leftwinvs(sub,:) = ep.icawinv(:,comps{sub}(1)) ; 
rightwinvs(sub,:) = ep.icawinv(:,comps{sub}(2)) ; 
end

subplot(2,2,1) ; imagesc(uniques,freqs,squeeze(mean(lefts,1)),[-5,5]) ; axis xy ; xlabel('polar angle(deg)') ; ylabel('frequency(hz)') ; 
subplot(2,2,2) ; topoplot(mean(leftwinvs,1),ep.chanlocs) ;

subplot(2,2,3) ; imagesc(uniques,freqs,squeeze(mean(rights,1)),[-5,5]) ; axis xy ; xlabel('polar angle(deg)') ; ylabel('frequency(hz)') ; 
subplot(2,2,4) ; topoplot(mean(rightwinvs,1),ep.chanlocs) ;

subplot(2,1,1) ; 
plot(smooth(squeeze(mean(mean(lefts(:,freqs>8 & freqs<15,:),2),1))),'LineWidth',3) ; hold on ; plot(smooth(squeeze(mean(mean(rights(:,freqs>8 & freqs<15,:),2),1))),'r','LineWidth',3) ; hline(0,'k') ; ylim([-6,2]) ; 
legend({'left component','right component'}) ; title('8-15Hz') ; set(gca,'XTick',1:5:117,'XTickLabel',round((1:5:117)*3)) ; xlabel('polar angle(deg)') ; ylabel('power(db)') ; 
subplot(2,1,2) ; 
plot(smooth(squeeze(mean(mean(lefts(:,freqs>40 & freqs<60,:),2),1))),'LineWidth',3) ; hold on ; plot(smooth(squeeze(mean(mean(rights(:,freqs>40 & freqs<60,:),2),1))),'r','LineWidth',3) ; hline(0,'k') ; ylim([-2,2]) ;
title('40-60Hz') ; set(gca,'XTick',1:5:length(uniques),'XTickLabel',round(uniques(1:5:end))) ; title('8-15Hz') ; set(gca,'XTick',1:5:117,'XTickLabel',round((1:5:117)*3)) ; 



%for i=1:64 ; subplot(5,13,i) ; imagesc(imfilter(squeeze(anglepower(i,:,:)),fspecial('gaussian',3,1)),[-3,3]) ; title(i) ; end
%for i=1:64 ; smoothpow(i,:,:) =
%imfilter(squeeze(anglepower(i,:,:)),fspecial('gaussian',3,1)) ; end
%figure,plot(squeeze(mean(smoothpow(25,freqs>10 & freqs<15,:),2))) ; hold on ; plot(squeeze(mean(smoothpow(25,freqs>40 & freqs<60,:),2)),'r') ; hline(0,'k')
%plot(squeeze(mean(smoothpow(24,freqs>10 & freqs<15,:),2)),squeeze(mean(smoothpow(24,freqs>40 & freqs<60,:),2)),'o') ;

