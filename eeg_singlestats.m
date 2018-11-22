clear all ; close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','jeremie','julie','katrine','lisa','marc','marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','vincent'} ;
allsis = load('c:/shared/allicas/allsis.mat') ; allsis = allsis.allsis ; 
for s=1:length(subs) ; 
    cd(['c:/shared/allicas/',subs{s}]) ; ls 
    EEG = pop_loadset('ica_notch85.set') ; 
    eegs{1} = pop_epoch(EEG,{'S 11'},[-.85,2.85]) ; 
    eegs{2} = pop_epoch(EEG,{'S 12'},[-.85,2.85]) ; 
    eegs{3} = pop_epoch(EEG,{'S 13'},[-.85,2.85]) ; 
    eegs{4} = pop_epoch(EEG,{'S 14'},[-.85,2.85]) ; 
    eegs{5} = pop_epoch(EEG,{'S 15'},[-.85,2.85]) ; 
    eegs{6} = pop_epoch(EEG,{'S 16'},[-.85,2.85]) ; 
    clear ersp ; 
    for eeg=1:length(eegs) 
        goodcomps = eegs{eeg}.icaact(allsis(s,1:3),:,:) ; 
        for i=1:size(goodcomps,1) ; 
            for j=1:135
                [ersp(eeg,i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(goodcomps(i,:,j)),eegs{1}.pnts,[eegs{1}.xmin,eegs{1}.xmax],eegs{1}.srate,0,...
                    'plotersp','off','plotitc','off','baseline',NaN,'freqs',[1,120]) ; 
            end
        end
    end
    baset = find(times<0) ; normersp = zeros(size(ersp)) ;
    for i=1:size(ersp,1)
        for j=1:size(ersp,2)
            for k=1:size(ersp,3)
                baseijk = repmat(squeeze(mean(ersp(i,j,k,:,baset),5)),[1,size(ersp,5)]) ;
                normersp(i,j,k,:,:) = squeeze(ersp(i,j,k,:,:))-baseijk ; 
            end
        end
    end   
    
    mtersp = squeeze(mean(normersp(:,:,:,:,times>0 & times<2),5)) ; 
    clear sigs ; 
    for i=1:6
        for j=1:3
            for k=1:120
                sigs(i,j,k) = anova1(squeeze(mtersp([1,i],j,:,k))',[],'off') ; 
            end
        end
    end
    allsigs(s,:,:,:) = sigs ; 
    allersp(s,:,:,:) = squeeze(mean(mean(normersp(:,:,:,:,times>0 &times<2),3),5)) ; 
end

%for i=1:6 ; subplot(2,3,i) ; imagesc(squeeze(mean(mean(allnormersp(5,i,:,:,:,:),3),4)),[-3,3]) ; end


%{
plot((median(freqdiffs)),'LineWidth',3) ; set(gca,'XTick',1:10:length(freqs),'XTickLabel',freqs(1:10:end)) ; title('median p-values, probability of mean(60%rnd) != mean(unperturbed)') ; 
subplot(2,2,2) ; 
bar(sort(median(freqdiffs(:,freqs<7& freqs>0),2))) ;
hline(0.05,'r','p=0.05') ; title('THETA 14/22 significantly less for randomization') ; 
subplot(2,2,1) ; 
bar(sort(median(freqdiffs(:,freqs<75& freqs>65),2))) ;
hline(0.05,'r','p=0.05') ; title('GAMMA 14/22 significantly less for randomization') ; 

ecount = 1 ; 
for i=1:size(mersp,2) ;
    tmersp(:,:,:,i) = squeeze(mersp(:,i,:,:)) ; 
end

save_nii(make_nii(tmersp),'tmersp.nii.gz') ;

%}

mtersp = allersp ; 
clear allinds ; 
for i=1:22;
    for j=1:3
        maxind = find(squeeze(mtersp(i,1,j,35:end))==max(squeeze(mtersp(i,1,j,35:end)))) + 35 ; 
        allinds(i,j,:) = (maxind-3):(maxind+3) ; 
    end
end

for i=1:22 ;
    for j=1:3
        sigdiffs(i,j) = median(squeeze(allsigs(i,2,j,squeeze(allinds(i,j,:))))) ; 
    end
end
bar(sort(min(sigdiffs,[],2))) ; hline(0.05,'r',{'p=0.05'}) ; xlabel('subjects') ; 
title('significant differences between contrast and unperturbed') ; 







