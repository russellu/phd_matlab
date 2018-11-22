clear all ; close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','jeremie','julie','katrine','lisa','marc','marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','vincent'} ;
allsis = load('c:/shared/allicas/allsis.mat') ; allsis = allsis.allsis ; 
stims = [2,3,1,5,6] ; 
for s=1:length(subs) ; 
    cd(['c:/shared/allres/',subs{s}]) ; ls 
    EEG = pop_loadset('ica_notch85.set') ; 
    %EEG = pop_subcomp(EEG,allsis(s,4:end)) ; 
    %EEG.data = abs(eegfiltfft(EEG.data,EEG.srate,50,80)) ; 
    ag = load('allgoods.mat') ; ag = ag.allgoods ; 
    eegs{1} = pop_epoch(EEG,{'S 11'},[-.85,2.85]) ; 
    eegs{2} = pop_epoch(EEG,{'S 12'},[-.85,2.85]) ; 
    eegs{3} = pop_epoch(EEG,{'S 13'},[-.85,2.85]) ; 
    eegs{4} = pop_epoch(EEG,{'S 14'},[-.85,2.85]) ; 
    eegs{5} = pop_epoch(EEG,{'S 15'},[-.85,2.85]) ; 
    eegs{6} = pop_epoch(EEG,{'S 16'},[-.85,2.85]) ;
    
    
    
        for e=1:length(eegs)
        for c=1:6 ; 
            for t=1:135
            [allersp(s,e,c,t,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(eegs{e}.icaact(allsis(s,c),:,t)),eegs{e}.pnts,[eegs{e}.xmin,eegs{e}.xmax],eegs{e}.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'baseline',NaN,'winsize',60) ; 
            end
        end   
        end
    %mersp = squeeze(mean(allersp(:,:,1:2,:,:),3)) ; 
    
    
    
    %{
    for e=1:length(eegs)
        for c=1:4 ; 
            [allersp(s,e,c,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(eegs{e}.icaact(allsis(s,c),:,ag{e})),eegs{e}.pnts,[eegs{e}.xmin,eegs{e}.xmax],eegs{e}.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'baseline',0,'winsize',60) ; 
            
        end   
    end
    mersp = squeeze(mean(allersp(:,:,1:2,:,:),3)) ; 
%}
    %{
    baset = 1:round(0.85*EEG.srate) ;
    taskt = round(0.85*EEG.srate):round(2.85*EEG.srate) ; 
    for i=1:length(eegs)
       tdiffs(i,:,:) = squeeze(mean(eegs{i}.icaact(:,taskt,1:135),2)-mean(eegs{i}.icaact(:,baset,1:135),2)) ;         
    end
    mdiffs = squeeze(mean(tdiffs,3)) ; 
    figure
    for i=1:6 ; subplot(2,3,i) ; 
        topoplot(double(mdiffs(i,:)),EEG.chanlocs,'maplimits',[-0.5,0.5]) ; 
    end
    alltdiffs(s,:,:,:) = tdiffs ; 
    disp(subs{s}) ; 
    %}
end
%{
stims = [2,3,1,5,6] ; 
for t=1:length(stims) ;
    msdiffs = squeeze(mean(mean(mean(alltdiffs(:,stims(t),:,:),1),2),4)) ; 
    subplot(1,5,t) ; 
    topoplot(double(msdiffs),EEG.chanlocs,'maplimits',[-0.3,0.3]) ; 
end
%}

dfs = find(freqs>6) ; 
t = [42,101,160] ;
f = [20,40,60,80,100] ; 
finds = [8,17,27,37,47]+3 ; 
for i=1:length(stims) ; 
    subplot(1,5,i) ; 
    imagesc(flipud(squeeze(mean(mersp(:,stims(i),dfs,:),1))),[-1.6,1.6]) ; colormap hot ; 
    set(gca,'XTick',t,'XTickLabel',[0,1,2],'YTick',finds,'YTickLabel',fliplr(f)) ; 
    if i==1
         xlabel('time(s)')  ; ylabel('frequency(hz)') ; 
    end
end

mtersp = squeeze(mean(mersp(:,:,:,times>0 & times<2),4)) ; 
errorbar(squeeze(mean(mtersp(:,[1,4],:),1))',squeeze(std(mtersp(:,[1,4],:),0,1))'./sqrt(22)) ; hline(0,'k') ; 

for i=1:60 
    for j=1:200
        [p,atab,stats] = anova1(squeeze(mersp(:,[1,4],i,j)),[],'off') ; 
        rsig(i,j) = atab{2,5} ; 
        
    end
end
imagesc(flipud(rsig)) ; colormap hot







mersp = squeeze(mean(mean(allersp,1),3)) ;
mtersp = squeeze(mean(mersp(:,:,times>0 & times<2),3)) ;
clear sv si
for i=1:6
   [sv(i,:),si] = sort(mtersp(i,freqs>35),'descend') ;  
   
end






msubs = squeeze(mean(mean(allersp(:,:,1:2,:,times>0 & times<2),5),3)) ; 
for i=1:6 ; subplot(2,3,i) ; imagesc(squeeze(msubs(:,i,:))) ; end
dfs = find(freqs>6) ; 
f = [20,40,60,80,100] ; 
finds = [8,17,27,37,47]+3 ; 
[sv,si] = sort(squeeze(mean(mean(msubs(:,:,freqs>35),2),3)),'descend') ; 
for i=1:5 ; 
    subplot(2,5,i) ; imagesc(fliplr(squeeze(msubs(si,stims(i),dfs)))',[-3,3]) ;
     set(gca,'YTick',finds,'YTickLabel',fliplr(f)) ; 
    if i==1 ; xlabel('subjects (n=22)') ; ylabel('frequency(hz)') ; end

end
subplot(2,5,10) ; imagesc([-3,3]);  colorbar ; 







bar(squeeze(mean(msubs(:,1,freqs>60 & freqs<80),3) - mean(msubs(:,2,freqs>60 & freqs<80),3)))





clear sdiffs ; for i=1:22 ; sdiffs(i,:,:) = squeeze(mean(alltdiffs(i,:,allsis(i,1:3),:),3)) ; end



mersp = squeeze(mean(allersp,3)) ; 

for i=1:22 ; disp(i)  ;
   for j=1:60
      for k=1:200         
           [h(i,j,k),p(i,j,k),ci,stats] = ttest(squeeze(mersp(i,1,:,j,k))',squeeze(mersp(i,6,:,j,k))') ; 
           subtvals(i,j,k) = stats.tstat ; 
      end
   end
end
for i=1:5 ; subplot(4,6,i) ; imagesc(squeeze(subtvals(i,:,:))) ; end
