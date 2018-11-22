clear all ; close all 
cd C:\shared\allres
subs = {'alex','alexandra3','audrey','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','leila','lisa','marc',...
    'marie','mathieu','maxime','menglu','mingham','olga','patricia','po','russell','suhan2','sunachakan','tah','tegan2','vincent'} ; 

for sub=1:6%length(subs)
    cd(['c:/shared/allres/',subs{sub}]) ; 
    trigs = dir('trig_*set') ; 
    trig = pop_loadset(trigs(1).name) ; 
    freqs = 1:2:100 ; 
    freqstims = zeros(length(freqs),size(trig.data,1),size(trig.data,2),size(trig.data,3)) ; 
    for f=1:length(freqs)
        for e=1:64
           freqstims(f,e,:,:) = (abs(eegfiltfft(squeeze(trig.data(e,:,:))',trig.srate,freqs(f)-1,freqs(f)+1)')) ; 
        end
    end    
    clear basestims ; 
    %basestims = freqstims - repmat(mean(freqstims(:,:,1:200,:),3),[1,1,1024,1]) ; 
    mgamma = squeeze(mean(mean(freqstims(:,:,300:700,:),3),4)) ; 
    allgamma(sub,:,:) = (mgamma) ; 
    figure,subplot(1,2,1) ; topoplot(mean(mgamma(32:37,:)),trig.chanlocs) ; subplot(1,2,2) ; imagesc(log(mgamma)) ; 
end

for i=1:6 ; subplot(2,3,i) ; topoplot(squeeze(mean(allgamma(i,5:12,:))),trig.chanlocs,'maplimits',[-2,2]) ; end


