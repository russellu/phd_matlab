cd c:/shared/allres
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','jeremie','julie','katrine','lisa'...
        ,'marc','marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','vincent'} ;

for s=1:length(subs) ;
   cd(['c:/shared/allres/',subs{s}]) ;  
   EEG = pop_loadset('ica_notch85.set') ; ag = load('allgoods') ; ag = ag.allgoods ; 
   trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'} ;
   for t=1:length(trigs)
        epst = pop_epoch(EEG,{trigs{t}},[-.85,2.85]) ; 
        for ch=1:size(epst.data,1)
           [erspi(s,t,ch,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epst.data(ch,:,ag{t})),epst.pnts,[epst.xmin,epst.xmax],epst.srate,0,...
               'plotersp','off','plotitc','off','winsize',64,'baseline',0,'freqs',[1,120],'nfreqs',60) ; 
        end
   end
end

labs = parse_freesurfer ; 
meaneeg = squeeze(mean(mean(erspi(:,:,:,:,times>0&times<2),2),5)) ; 
%meanf = squeeze(mean(mean(allstimepochs(:,:,:,5:6),3),4)) ; meaneeg = meanf ; 
clear corrmat ; 
for i=1:size(labs,2)
    labi = cell2mat({labs{:,i,2}}) ;
    for j=1:size(labi,1) ; 
        for k=1:size(meaneeg,2)
            for el=1:size(meaneeg,3)
                corrmat{i}(j,k,el) = corr2(meaneeg(:,k,el),labi(j,:)') ;
            end
        end
    end
end
es = {EEG.chanlocs.labels} ; 
for e=1:64 ; figure,
for i=1:19 ; subplot(4,5,i) ; imagesc(squeeze(corrmat{i}(:,e,:)),[-1,1]) ; end
suptitle(es{e}) ; 
end
imagesc(corrmat{9},[-1,1]) ; set(gca,'YTick',1:size(corrmat{1},1),'YTickLabel',labs{1,9,1}) ; 
param =9; hzind = 5; paramrow = 13 ; 
labi = cell2mat({labs{:,param,2}}) ;
figure,plot(squeeze(meaneeg(:,hzind)),labi(paramrow,:),'o') ; title(num2str(corr2(squeeze(meaneeg(:,hzind))',labi(paramrow,:)))) ; 
figure,imagesc(corrmat{param},[-1,1]) ; set(gca,'YTick',1:size(corrmat{param},1),'YTickLabel',labs{1,param,1}) ; 






for i=1:size(corrmat,2) ; figure ; 
    for j=1:size(corrmat{i},1) ;
        subplot(5,10,j) ; topoplot(squeeze(mean(corrmat{i}(j,:,freqs>30 & freqs<50),3)),EEG.chanlocs,'maplimits',[-1,1]) ; 
        title(labs{1,i,1}(j)) ; 
    end ; 
end









