clear all ; close all ; 
responders = {'alexandra3','fabio','gab','gabriella','genevieve','gina','jeremie','julie','katrine','marie','maxime','mingham','po','russell','suhan2','tegan2'} ;
stims = {'S 11','S 12','S 13','S 14','S 15','S 16'} ;
% get the top components (hand picked)
for r=16%:length(responders)
   cd(['c:/shared/allres/',responders{r}]) ; ls   
   EEG = pop_loadset('ica_notch85.set') ; 
   stimset = pop_epoch(EEG,{stims{[1]}},[-1,3]) ; 
   for e=1:64
   for i=1:size(stimset.data,3)
      [ersp(e,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(stimset.icaact(e,:,i)),stimset.pnts,[stimset.xmin,stimset.xmax],stimset.srate,0,...
          'plotersp','off','plotitc','off','nfreqs',120,'freqs',[1,120],'winsize',40,'baseline',NaN) ; 
   end
   end
    
end

for i=1:size(ersp,1)
    for j=1:size(ersp,2)
        for k=1:size(ersp,3)
            baseijk = repmat(squeeze(mean(ersp(i,j,k,times<0),4)),[1,200]) ; 
            bersp(i,j,k,:) = squeeze(ersp(i,j,k,:))' - baseijk ; 
        end
    end
end



bursts = (squeeze(mean(bersp(10,:,freqs>60&freqs<80,times>0.1 & times<2.1),3))) ;
fbursts = (squeeze(mean(bersp(10,:,freqs>60&freqs<80,:),3))) ;

cmat = corr(bursts') ; for i=1:size(cmat,1) ; for j=1:size(cmat,2) ; if i<=j cmat(i,j) = 0 ; end ; end ; end
fcmat = corr(bursts') ; 
corrvals = cmat(cmat~=0) ; posinds = find(cmat>0.4) ; 
clear w h
for i=1:length(posinds)
    [w(i),h(i)] = ind2sub(size(cmat),posinds(i)) ; 
    
end
for i=1:50 ; 
    subplot(5,10,i) ;
    plot(bursts(w(i),:)-mean(bursts(w(i),:)),'LineWidth',1) ; hold on ; 
    plot(bursts(h(i),:)-mean(bursts(h(i),:)),'r','LineWidth',1) ; xlim([0,100]) ; 
    title(num2str(cmat(w(i),h(i)))) ; 
end

for s=1:50 ; figure,
subplot(1,2,1) ; imagesc(squeeze(bersp(15,w(s),:,:)),[-15,15]) ; subplot(1,2,2) ; imagesc(squeeze(bersp(15,h(s),:,:)),[-15,15])
end



for i=1:size(fcmat,1)
   inds{i} = find(fcmat(i,:)>0.4) ;  
    
    
end

for i=1:length(inds)
    subplot(10,14,i) ; imagesc(squeeze(mean(bersp(10,inds{i},:,:),2)),[-10,10]) ; 
    
end

